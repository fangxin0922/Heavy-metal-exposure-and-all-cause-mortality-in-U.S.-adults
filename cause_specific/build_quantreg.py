#!/usr/bin/env python3
"""Build the study's isolated quantreg 6.1 sparse-event memory repair."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,tarfile,urllib.request,urllib.error

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'cause_specific'
ARCHIVE_SHA='f42292c5ab25a15f39295b93391deafef192fe09eefde563399a64eba7e0169a'
ORIGINAL_SHA='b45374d4e186df6ce72191cc1e4862a293827dda2a30aa55d3b0e1a57d653fbd'
REPAIRED_SHA='c934cf702f62d75d34d0d475221a8bebf47fffcff83f0b90a14c29c6823121a1'

def download_source(destination):
    # CRAN moves a version to Archive only after it is superseded.
    urls=['https://cran.r-project.org/src/contrib/quantreg_6.1.tar.gz',
          'https://cran.r-project.org/src/contrib/Archive/quantreg/quantreg_6.1.tar.gz']
    for url in urls:
        try:
            urllib.request.urlretrieve(url,destination)
            return url
        except urllib.error.HTTPError as error:
            if error.code not in (404,410):raise
    raise RuntimeError('Pinned quantreg source unavailable from CRAN; supply --archive.')

def repair(source):
    if hashlib.sha256(source).hexdigest()!=ORIGINAL_SHA:
        raise ValueError('Unexpected quantreg native source')
    lines=source.decode().splitlines(keepends=True)
    k=next(i for i,line in enumerate(lines) if line.rstrip()=='      do23002 i = 1,n')
    assert lines[k+1].rstrip()=='      u(i) = one' and lines[k+2].rstrip()=='      wn(i,1) = half'
    result=''.join(lines[:k]+['      do23004 i = 1,n1\n','      wn(i,1) = half\n','23004 continue\n']+lines[k:k+2]+lines[k+3:]).encode()
    assert hashlib.sha256(result).hexdigest()==REPAIRED_SHA
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive',type=Path,help='Existing unmodified quantreg_6.1.tar.gz')
    args=parser.parse_args()
    build=OUT/'build';library=OUT/'library'
    if (library/'quantreg').exists() or (build/'quantreg').exists():
        raise RuntimeError('Use an empty cause_specific build/library directory.')
    build.mkdir(exist_ok=True);library.mkdir(exist_ok=True)
    archive=build/'quantreg_6.1.tar.gz'
    if args.archive: shutil.copyfile(args.archive,archive)
    else: download_source(archive)
    if hashlib.sha256(archive.read_bytes()).hexdigest()!=ARCHIVE_SHA:
        raise RuntimeError('The quantreg source archive does not match the manuscript version.')
    with tarfile.open(archive) as tar: tar.extractall(build,filter='data')
    native=build/'quantreg/src/crqfnb.f'
    native.write_bytes(repair(native.read_bytes()))
    env=os.environ.copy()
    if env.get('HM_R_LIBRARY'):env['R_LIBS_USER']=env['HM_R_LIBRARY']
    with (build/'install.log').open('w') as log:
        subprocess.run(['R','CMD','INSTALL','--preclean','--no-multiarch','--no-docs','--no-html',
                        '-l',str(library),str(build/'quantreg')],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    subprocess.run(['Rscript','cause_specific/R/record_build.R'],cwd=ROOT,env=env,check=True)
    (build/'source_provenance.json').write_text(json.dumps(dict(archive_sha256=ARCHIVE_SHA,
        original_native_sha256=ORIGINAL_SHA,repaired_native_sha256=REPAIRED_SHA,
        source_change='Initialise wn over n1, separately from u over n; no estimator or grid changes'),indent=2)+'\n')
    print('Isolated quantreg repair installed:',library)

if __name__=='__main__':main()
