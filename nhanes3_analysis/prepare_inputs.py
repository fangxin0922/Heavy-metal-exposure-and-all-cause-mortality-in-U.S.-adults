#!/usr/bin/env python3
"""Download verified NHANES III public releases and extract documented ASCII fields."""
from pathlib import Path
import argparse,csv,hashlib,json,os,re,shutil,urllib.request
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source-dir',type=Path,help='Existing original raw files; otherwise download from CDC')
args=parser.parse_args()
root=Path(os.environ.get('NH3_OUTPUT_DIR',str(Path(__file__).resolve().parent)))
root.mkdir(parents=True,exist_ok=True)
raw=root/'raw';raw.mkdir(exist_ok=True)
sources=json.loads((Path(__file__).resolve().parent/'sources.json').read_text())
for item in sources:
    target=raw/item['file']
    if not target.exists():
        temporary=target.with_suffix(target.suffix+'.partial')
        if args.source_dir: shutil.copyfile(args.source_dir/item['file'],temporary)
        else:
            request=urllib.request.Request(item['url'],headers={'User-Agent':'NHANES-reproduction/1.0'})
            with urllib.request.urlopen(request,timeout=300) as response,temporary.open('wb') as out:
                shutil.copyfileobj(response,out)
        if hashlib.sha256(temporary.read_bytes()).hexdigest()!=item['sha256']:
            raise RuntimeError('Public input differs from the documented release: '+item['file'])
        temporary.replace(target)
    if hashlib.sha256(target.read_bytes()).hexdigest()!=item['sha256']:
        raise RuntimeError('Input checksum mismatch: '+str(target))
(raw/'verified_sources.json').write_text(json.dumps(sources,indent=2)+'\n')
selected={
'exam': 'SEQN DMPSTAT HSSEX HSAGEIR DMARETHN DMARACER DMAETHNR DMPPIR SDPPHASE SDPSTRA6 SDPPSU6 WTPFEX6 BMPBMI MAPE1 MAPE2 MAPE3S MAPE4 MAPF12R'.split(),
'lab': 'SEQN PBP UDP URP CEP TCP'.split(),
'adult': 'SEQN HFA8R HFA12 HAR1 HAR2 HAR3 HAR4S HAR9 HAR10 HAR11R HAR12S HAT1S HAT2 HAT4 HAT6 HAT8 HAT10 HAT12 HAT14 HAT16 HAT18'.split(),
}
(root/'data').mkdir(exist_ok=True)
if (root/'data/exam_selected.csv').exists(): raise RuntimeError('Use an empty NHANES III data directory.')
audit=[]
for name,wanted in selected.items():
    sas=(root/'raw'/f'{name}.sas').read_text()
    block=re.search(r'\bINPUT\s+(.*?);',sas,re.S).group(1)
    positions={m[1]:(int(m[2]),int(m[3] or m[2])) for m in re.finditer(r'^\s*(\w+)\s+(\d+)(?:-(\d+))?\s*$',block,re.M)}
    assert set(wanted)<=positions.keys(),set(wanted)-positions.keys()
    # Official ADULTMRG drops HAX* and adds ADULTX: none of these selected fields is replaced.
    if name=='adult':
        update=(root/'raw'/'adultx.sas').read_text()
        ub=re.search(r'\bINPUT\s+(.*?);',update,re.S).group(1)
        uvars=set(re.findall(r'^\s*(\w+)\s+\d',ub,re.M))
        assert not (set(wanted)-{'SEQN'}) & uvars
        assert not any(v.startswith('HAX') for v in wanted)
    ids=set();rows=0
    with (root/'data'/f'{name}_selected.csv').open('w') as out:
        writer=csv.writer(out);writer.writerow(wanted)
        for line in (root/'raw'/f'{name}.dat').open():
            vals=[line[positions[v][0]-1:positions[v][1]].strip() for v in wanted]
            assert vals[0] and vals[0] not in ids
            ids.add(vals[0]);writer.writerow(vals);rows+=1
    for v in wanted:audit.append(dict(source=name,variable=v,start=positions[v][0],end=positions[v][1],rows=rows))
    print(name,rows,flush=True)
positions={'SEQN':(1,6),'ELIGSTAT':(15,15),'MORTSTAT':(16,16),'PERMTH_EXM':(46,48)}
with (root/'data'/'mortality_selected.csv').open('w') as out:
    w=csv.writer(out);w.writerow(positions);ids=set()
    for line in (root/'raw'/'mortality.dat').open():
        vals=[line[a-1:b].strip() for a,b in positions.values()]
        assert vals[0] not in ids;ids.add(vals[0]);w.writerow(vals)
    print('mortality',len(ids),flush=True)
(root/'data'/'extraction_dictionary.json').write_text(json.dumps(audit,indent=2)+'\n')
