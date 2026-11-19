#!/usr/bin/env python3
"""Run cause-censored CQR using the completed main primary models."""
from pathlib import Path
import argparse,fcntl,os,sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from pipeline_runner import ROOT,execute,parallel

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--points-only',action='store_true')
    parser.add_argument('--fit-only',action='store_true')
    args=parser.parse_args()
    if not 1<=args.workers<=24:parser.error('workers must be 1–24')
    out=Path(os.environ.get('CAUSE_OUTPUT_DIR',str(ROOT/'cause_specific')))
    if not out.is_absolute():out=ROOT/out
    out.mkdir(parents=True,exist_ok=True)
    with (out/'RUN.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if not args.fit_only:execute(['Rscript','--vanilla','cause_specific/R/prepare_models.R'],out/'logs/prepare.log')
        scenarios=['primary_Pb','primary_Hg','egfr60_Cd','egfr60_As']
        parallel([(['Rscript','--vanilla','cause_specific/R/fit_models.R',ep,sc,str(i),'0' if args.points_only else '400'],out/f'logs/{ep}_{sc}_{i}.log') for ep in ['cancer','cardiovascular'] for sc in scenarios for i in range(1,6)],args.workers)
        execute(['Rscript','--vanilla','scripts/summarize_extension.R','cause'],out/'logs/summarize.log')

if __name__=='__main__':main()
