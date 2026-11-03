"""Fit biomarker interaction models with bounded, single-threaded R workers."""
import argparse
import concurrent.futures
import csv
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("original", "replicates"))
    parser.add_argument("workers", type=int, nargs="?", default=24)
    args = parser.parse_args()
    if not 1 <= args.workers <= 24:
        parser.error("workers must be between 1 and 24")

    root = Path(__file__).resolve().parent
    output = Path(os.environ.get("EMM_OUTPUT_DIR", "biomarker_emm/output"))
    if not output.is_absolute():
        output = root.parent / output
    with (output / "aggregate/model_preparation.csv").open() as handle:
        pairs = [row["pair"] for row in csv.DictReader(handle)]
    assert len(pairs) == len(set(pairs)) == 12
    blocks = ([(0, 0)] if args.stage == "original" else
              [(i, min(i + 19, 400)) for i in range(1, 401, 20)])
    jobs = [(pair, first, last) for pair in pairs for first, last in blocks]
    (output / "logs").mkdir(parents=True, exist_ok=True)
    (output / "aggregate").mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

    def run(job):
        pair, first, last = job
        log = output / "logs" / f"{pair}_{first:03d}_{last:03d}.log"
        command = ["Rscript", str(root / "R/fit_models.R"), pair, str(first), str(last)]
        with log.open("a") as stream:
            result = subprocess.run(command, cwd=root.parent, env=env,
                                    stdout=stream, stderr=subprocess.STDOUT)
        return dict(pair=pair, first=first, last=last,
                    exit_code=result.returncode, log=str(log))

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(run, jobs))
    with (output / f"aggregate/{args.stage}_jobs.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    failed = sum(row["exit_code"] != 0 for row in results)
    print(f"{args.stage}: {len(results) - failed}/{len(results)} jobs finished successfully", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
