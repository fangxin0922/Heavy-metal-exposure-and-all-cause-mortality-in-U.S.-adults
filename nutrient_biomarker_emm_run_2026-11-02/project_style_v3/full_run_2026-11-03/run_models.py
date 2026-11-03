"""Run bounded, restartable R jobs; participant inputs stay on the server."""
import concurrent.futures
import csv
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path("nutrient_biomarker_emm_run_2026-11-02/project_style_v3")
RUN_DIR = ROOT / "full_run_2026-11-03"
STAGE = sys.argv[1]
WORKERS = int(sys.argv[2]) if len(sys.argv) > 2 else 24
assert STAGE in {"original", "pilot", "replicates"}
assert 1 <= WORKERS <= 24
with (ROOT / "aggregate/model_preparation.csv").open() as handle:
    pairs = [row["pair"] for row in csv.DictReader(handle)]
assert len(pairs) == len(set(pairs)) == 12
blocks = [(0, 0)] if STAGE == "original" else ([(i, min(i + 3, 20)) for i in range(1, 21, 4)] if STAGE == "pilot" else [(i, min(i + 19, 400)) for i in range(1, 401, 20)])
jobs = [(pair, first, last) for pair in pairs for first, last in blocks]
(RUN_DIR / "logs").mkdir(exist_ok=True)
env = os.environ.copy()
env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

def run(job):
    pair, first, last = job
    log = RUN_DIR / "logs" / f"{pair}_{first:03d}_{last:03d}.log"
    command = ["Rscript", str(ROOT / "R/03_fit_cqr.R"), pair, str(first), str(last)]
    with log.open("a") as output:
        result = subprocess.run(command, env=env, stdout=output, stderr=subprocess.STDOUT)
    return dict(pair=pair, first=first, last=last, exit_code=result.returncode, log=str(log))

with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as pool:
    results = list(pool.map(run, jobs))
with (RUN_DIR / f"aggregate/{STAGE}_jobs.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)
failed = sum(row["exit_code"] != 0 for row in results)
print(f"{STAGE}: {len(results) - failed}/{len(results)} jobs finished successfully", flush=True)
sys.exit(1 if failed else 0)
