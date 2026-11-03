#!/usr/bin/env python3
"""Schedule the R analysis and track each job."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import argparse
import csv
import datetime
import fcntl
import hashlib
import json
import os
import subprocess
import threading
import time

ROOT = Path(__file__).resolve().parent.parent
RUN = Path(__file__).resolve().parent
REL = RUN.relative_to(ROOT)
PRIMARY = ['primary_Pb', 'primary_Hg', 'egfr60_Cd', 'egfr60_As']
ALL = PRIMARY + ['primary_Cd', 'primary_As']
ALL += [f'{panel}_{metal}' for panel in ['no_bmi', 'age55plus', 'portnoy']
        for metal in ['Pb', 'Hg', 'Cd', 'As']]
ALL += ['joint_Pb', 'joint_Cd', 'epa_dha_Hg', 'rcs_Pb', 'rcs_Cd']
ENV = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           MKL_NUM_THREADS='1', BLIS_NUM_THREADS='1')

def now():
    return datetime.datetime.now().astimezone().isoformat()

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=16)
    parser.add_argument('--prepare-workers', type=int, default=4)
    parser.add_argument('--scope', choices=['all', 'primary'], default='all')
    parser.add_argument('--points-only', action='store_true', help='Explicit diagnostic mode, no intervals')
    args = parser.parse_args()
    assert 1 <= args.workers <= 24 and 1 <= args.prepare_workers <= 8
    scenarios = ALL if args.scope == 'all' else PRIMARY
    assert len(ALL) == 23 and len(set(ALL)) == 23
    (RUN / 'logs').mkdir(parents=True, exist_ok=True)
    (RUN / 'audit').mkdir(parents=True, exist_ok=True)
    lock = open(RUN / 'RUN.lock', 'a+')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit('The fresh pipeline already has an active controller.')
    lock.seek(0)
    lock.truncate()
    lock.write(str(os.getpid()))
    lock.flush()
    source_files = list((RUN / 'R').glob('*.R')) + [Path(__file__).resolve()]
    sources = {str(p.relative_to(ROOT)): sha(p) for p in sorted(source_files)}
    manifest_path = RUN / 'RUN_MANIFEST.json'
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text())
        if previous['source_sha256'] != sources:
            raise SystemExit('Source changed since launch. Start a separately named run or review explicitly.')
    else:
        manifest_path.write_text(json.dumps(dict(created=now(), source_sha256=sources,
            reference_commit='da34a8bdeff9d665da38f7f07c853aaaf71bc5c3'), indent=2))
    with open(RUN / 'audit' / 'requested_jobs.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['scenario', 'imputation', 'resamples'])
        writer.writerows((scenario, mi, 0 if args.points_only else 400)
                        for scenario in scenarios for mi in range(1, 6))
    state = dict(started=now(), controller_pid=os.getpid(), status='running',
                 scope=args.scope, bootstrap_resamples=0 if args.points_only else 400,
                 preparation={}, fitting={}, stages={}, errors=[])
    mutex = threading.Lock()

    def update(section, name, value):
        with mutex:
            state[section][name] = value
            state['updated'] = now()
            temporary = RUN / 'RUN_STATUS.json.tmp'
            temporary.write_text(json.dumps(state, indent=2))
            temporary.replace(RUN / 'RUN_STATUS.json')

    def run_r(script, arguments, section, name):
        logfile = RUN / 'logs' / (name + '.log')
        with open(logfile, 'w') as log:
            process = subprocess.Popen(['Rscript', '--vanilla', str(REL / 'R' / script)]
                + [str(x) for x in arguments], cwd=ROOT, env=ENV, stdin=subprocess.DEVNULL,
                stdout=log, stderr=subprocess.STDOUT)
            update(section, name, dict(status='running', pid=process.pid, started=now(), log=str(logfile)))
            code = process.wait()
        update(section, name, dict(status='complete' if code == 0 else 'failed',
                                   returncode=code, finished=now(), log=str(logfile)))
        return code

    try:
        # Reuse inputs only after their hashes match.
        for number, script in enumerate(['exposure_records.R',
                'survival_records.R', 'adjustment_data.R',
                'cohort_builder.R'], 1):
            marker = RUN / 'audit' / f'stage{number:02d}.complete'
            if not marker.exists():
                if run_r(script, [], 'stages', f'input_{number:02d}'):
                    raise RuntimeError(f'Input stage failed: {script}')
                marker.write_text(now())
            else:
                update('stages', f'input_{number:02d}', dict(status='completed_in_this_run'))
        frozen = RUN / 'audit' / 'frozen_inputs.json'
        hashes = {p.name: sha(p) for p in [RUN / 'data' / 'analytic_fresh.rds', RUN / 'data' / 'full_designs.rds']}
        if frozen.exists():
            assert json.loads(frozen.read_text()) == hashes, 'Frozen fresh input changed'
        else:
            frozen.write_text(json.dumps(hashes, indent=2))

        if run_r('replay_baseline.R', [], 'stages', 'control_replay'):
            state['errors'].append('Control replay failed; inspect its job logs')
        run_r('compare_runs.R', [], 'stages', 'control_comparison')

        # Queue model fits after each scenario finishes imputation.
        fit_futures = []
        with ThreadPoolExecutor(max_workers=args.workers) as fitting:
            def prepare(scenario):
                code = run_r('impute_covariates.R', [scenario], 'preparation', scenario)
                if code:
                    return scenario, code
                for mi in range(1, 6):
                    out = RUN / 'fits' / 'fresh' / scenario / f'imp{mi:02d}'
                    # Validate saved results before reuse.
                    future = fitting.submit(run_r, 'fit_cqr.R',
                        [scenario, mi, 0 if args.points_only else 400, 'fresh'],
                        'fitting', f'{scenario}_{mi}')
                    with mutex:
                        fit_futures.append(future)
                return scenario, 0

            with ThreadPoolExecutor(max_workers=args.prepare_workers) as preparation:
                futures = [preparation.submit(prepare, scenario) for scenario in scenarios]
                for future in as_completed(futures):
                    scenario, code = future.result()
                    if code:
                        state['errors'].append(f'Preparation failed: {scenario}')
                    run_r('compare_runs.R', [], 'stages', 'progress_report')
            # Collect every job status, including failures.
            last_report = time.monotonic()
            for future in as_completed(fit_futures):
                if future.result():
                    state['errors'].append('A fitting process failed; see fitting status')
                if time.monotonic() - last_report > 60:
                    run_r('compare_runs.R', [], 'stages', 'progress_report')
                    last_report = time.monotonic()

        for script, name in [('compare_runs.R', 'final_comparison'),
                ('survey_cox.R', 'descriptive_cox'), ('draw_figures.R', 'figures')]:
            if run_r(script, [], 'stages', name):
                state['errors'].append(f'Stage failed: {script}')
        state['status'] = 'complete' if not state['errors'] else 'finished_with_errors'
    except Exception as error:
        state['status'] = 'failed'
        state['errors'].append(str(error))
        raise
    finally:
        state['finished'] = now()
        update('stages', 'controller', dict(status=state['status']))
        print(json.dumps(dict(status=state['status'], errors=state['errors']), indent=2), flush=True)

if __name__ == '__main__':
    main()
