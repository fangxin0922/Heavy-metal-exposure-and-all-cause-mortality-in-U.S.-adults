#!/usr/bin/env python3
"""Run the metal survival analysis from input preparation to figures."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import argparse
import fcntl
import os
import subprocess

ROOT = Path(__file__).resolve().parent.parent
RUN = Path(__file__).resolve().parent
REL = RUN.relative_to(ROOT)
PRIMARY = ['primary_Pb', 'primary_Hg', 'egfr60_Cd', 'egfr60_As']
ALL = PRIMARY + ['primary_Cd', 'primary_As']
ALL += [f'{panel}_{"egfr60_" if metal in ["Cd", "As"] else ""}{metal}'
        for panel in ['no_bmi', 'age55plus', 'portnoy']
        for metal in ['Pb', 'Hg', 'Cd', 'As']]
ALL += ['joint_egfr60_Pb', 'joint_egfr60_Cd', 'epa_dha_Hg', 'rcs_Pb', 'rcs_Cd']
PREPARATION_DEPENDENCIES = {
    f'{panel}_egfr60_{metal}': f'egfr60_{metal}'
    for panel in ['no_bmi', 'age55plus', 'portnoy'] for metal in ['Cd', 'As']}
PREPARATION_DEPENDENCIES.update({
    'joint_egfr60_Pb': 'egfr60_Cd', 'joint_egfr60_Cd': 'egfr60_Cd'})
ENV = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           MKL_NUM_THREADS='1', BLIS_NUM_THREADS='1')
OUTPUTS = ['data', 'models', 'fits', 'results', 'logs', 'diagnostics', 'figures']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=16)
    parser.add_argument('--prepare-workers', type=int, default=4)
    parser.add_argument('--scope', choices=['all', 'primary'], default='all')
    parser.add_argument('--points-only', action='store_true', help='Fit CQR points without intervals')
    args = parser.parse_args()
    if not (1 <= args.workers <= 24 and 1 <= args.prepare_workers <= 8):
        parser.error('Use 1–24 fit workers and 1–8 preparation workers.')
    scenarios = ALL if args.scope == 'all' else PRIMARY
    assert len(ALL) == 23 and len(set(ALL)) == 23
    with open(RUN / 'RUN.lock', 'a+') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit('An analysis is already running.')
        occupied = [name for name in OUTPUTS
                    if (RUN / name).exists() and any((RUN / name).iterdir())]
        if occupied:
            raise SystemExit('Use empty output folders for a new run: ' + ', '.join(occupied))
        for name in OUTPUTS:
            (RUN / name).mkdir(parents=True, exist_ok=True)

        def run_r(script, arguments, name):
            logfile = RUN / 'logs' / (name + '.log')
            print(f'Starting {name}', flush=True)
            with open(logfile, 'w') as log:
                result = subprocess.run(
                    ['Rscript', '--vanilla', str(REL / 'R' / script)] + list(map(str, arguments)),
                    cwd=ROOT, env=ENV, stdin=subprocess.DEVNULL,
                    stdout=log, stderr=subprocess.STDOUT, check=False)
            if result.returncode:
                raise RuntimeError(f'{name} failed; see {logfile}')
            print(f'Completed {name}', flush=True)

        for script in ['exposure_records', 'survival_records', 'adjustment_data', 'cohort_builder']:
            run_r(script + '.R', [], script)

        # Complete primary preparations before their dependent sensitivities;
        # primary fits may continue while the sensitivities are prepared.
        errors = []
        prepared = set()
        with ThreadPoolExecutor(max_workers=args.workers) as fitting:
            fit_jobs = []
            with ThreadPoolExecutor(max_workers=args.prepare_workers) as preparation:
                batches = [[scenario for scenario in scenarios if scenario in PRIMARY],
                           [scenario for scenario in scenarios if scenario not in PRIMARY]]
                for batch in batches:
                    preparation_jobs = {}
                    for scenario in batch:
                        dependency = PREPARATION_DEPENDENCIES.get(scenario)
                        if dependency is not None and dependency not in prepared:
                            errors.append(f'{scenario} skipped: required preparation '
                                          f'{dependency} failed.')
                            continue
                        future = preparation.submit(
                            run_r, 'impute_covariates.R', [scenario], scenario)
                        preparation_jobs[future] = scenario
                    for future in as_completed(preparation_jobs):
                        scenario = preparation_jobs[future]
                        try:
                            future.result()
                        except Exception as error:
                            errors.append(str(error))
                            continue
                        prepared.add(scenario)
                        for mi in range(1, 6):
                            fit_jobs.append(fitting.submit(
                                run_r, 'fit_cqr.R', [scenario, mi, 0 if args.points_only else 400],
                                f'{scenario}_{mi}'))
            for future in as_completed(fit_jobs):
                try:
                    future.result()
                except Exception as error:
                    errors.append(str(error))
        if errors:
            raise RuntimeError('\n'.join(errors))

        run_r('summarize_results.R', scenarios, 'summarize_results')
        run_r('survey_cox.R', [], 'survey_cox')
        run_r('draw_figures.R', [], 'draw_figures')
        print('Analysis complete. Tables are in results/ and figures are in figures/.', flush=True)


if __name__ == '__main__':
    main()
