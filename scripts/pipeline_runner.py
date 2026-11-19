"""Bounded R worker support for the two additional study pipelines."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import os,subprocess

ROOT=Path(__file__).resolve().parent.parent
ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')

def execute(command,log):
    log=Path(log);log.parent.mkdir(parents=True,exist_ok=True)
    print('Starting '+log.stem,flush=True)
    with log.open('w') as stream:
        result=subprocess.run(command,cwd=ROOT,env=ENV,stdin=subprocess.DEVNULL,stdout=stream,stderr=subprocess.STDOUT)
    if result.returncode: raise RuntimeError('Command failed; see '+str(log))
    print('Completed '+log.stem,flush=True)

def parallel(jobs,workers):
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(execute,command,log) for command,log in jobs]
        errors=[]
        for future in futures:
            try:future.result()
            except Exception as error:errors.append(str(error))
    if errors:raise RuntimeError('\n'.join(errors))
