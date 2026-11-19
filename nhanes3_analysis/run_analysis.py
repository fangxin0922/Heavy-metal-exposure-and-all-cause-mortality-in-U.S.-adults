#!/usr/bin/env python3
"""Run NHANES III public-input preparation, five imputations, and CQR."""
from pathlib import Path
import argparse,fcntl,os,sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from pipeline_runner import ROOT,execute,parallel

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir',type=Path)
    parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--points-only',action='store_true')
    parser.add_argument('--fit-only',action='store_true',help='Use already prepared NHANES III models')
    args=parser.parse_args()
    if not 1<=args.workers<=24:parser.error('workers must be 1–24')
    out=Path(os.environ.get('NH3_OUTPUT_DIR',str(ROOT/'nhanes3_analysis')))
    if not out.is_absolute():out=ROOT/out
    out.mkdir(parents=True,exist_ok=True)
    with (out/'RUN.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if not args.fit_only:
            command=[sys.executable,'nhanes3_analysis/prepare_inputs.py']
            if args.source_dir:command+=['--source-dir',str(args.source_dir.resolve())]
            execute(command,out/'logs/public_inputs.log')
            execute(['Rscript','--vanilla','nhanes3_analysis/R/build_cohorts.R'],out/'logs/build_cohorts.log')
            parallel([(['Rscript','--vanilla','nhanes3_analysis/R/impute_covariates.R',sc],out/f'logs/prepare_{sc}.log') for sc in ['lead','cadmium']],2)
        parallel([(['Rscript','--vanilla','nhanes3_analysis/R/fit_models.R',sc,str(i),'0' if args.points_only else '400'],out/f'logs/fit_{sc}_{i}.log') for sc in ['lead','cadmium'] for i in range(1,6)],args.workers)
        execute(['Rscript','--vanilla','scripts/summarize_extension.R','nhanes3'],out/'logs/summarize.log')

if __name__=='__main__':main()
