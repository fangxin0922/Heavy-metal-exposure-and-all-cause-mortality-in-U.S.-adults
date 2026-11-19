#!/usr/bin/env python3
"""Compare fresh outputs with manuscript reference tables; never replace fitted values."""
from pathlib import Path
import argparse,csv,json,math

ROOT=Path(__file__).resolve().parent.parent
def read(path):
    with path.open() as f:return list(csv.DictReader(f))
def number(value):
    try:return float(value)
    except (ValueError,TypeError):return float('nan')
def finite(value):return math.isfinite(number(value))
def output(path,rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

def compare(family,reference,fresh):
    keys=['pair','tau'] if family=='biomarker' else (['endpoint'] if family=='cause' else [])+['scenario','imputation','term','tau']
    def key(row):return tuple(str(float(row[k])) if k in ['tau','imputation'] else row[k] for k in keys)
    original={key(r):r for r in reference};new={key(r):r for r in fresh}
    if len(original)!=len(reference) or len(new)!=len(fresh):raise ValueError('Duplicate output cells: '+family)
    if original.keys()!=new.keys():raise ValueError('Missing or unexpected output cells: '+family)
    measures=[('ic_years','ic_lower_years','ic_upper_years'),('metal_low_years','metal_low_lower','metal_low_upper'),('metal_high_years','metal_high_lower','metal_high_upper')] if family=='biomarker' else [('estimate','lower','upper')]
    rows=[]
    for k,ref in original.items():
        rerun=new[k]
        for field in ['n','deaths']:
            if int(ref[field])!=int(rerun[field]):raise ValueError(f'{family} cohort/event mismatch: {k}: {field}')
        if family=='biomarker':
            for field in ['biomarker_low','biomarker_high']:
                if not math.isclose(number(ref[field]),number(rerun[field]),rel_tol=1e-10,abs_tol=1e-10):
                    raise ValueError(f'Biomarker reference level changed: {k}: {field}')
        for estimate,lower,upper in measures:
            a=number(ref[estimate]);b=number(rerun[estimate]);fa=math.isfinite(a);fb=math.isfinite(b)
            difference=b-a if fa and fb else None
            # A review flag, not a statistical equivalence margin. Preserve availability separately.
            threshold=max(.1,.2*abs(a)) if fa else None
            intervals=all(finite(row[field]) for row in [ref,rerun] for field in [lower,upper])
            inference_measure=family!='biomarker' or estimate=='ic_years'
            excludes=lambda row:number(row[lower])>0 or number(row[upper])<0
            rows.append(dict(family=family,cell='|'.join(k),measure=estimate,
                reference=a if fa else '',rerun=b if fb else '',difference_years=difference if difference is not None else '',
                same_availability=fa==fb,paired_finite=fa and fb,
                sign_changed=(a*b<0) if fa and fb else '',review_threshold_years=threshold if threshold is not None else '',
                magnitude_flag=abs(difference)>threshold if difference is not None else '',
                ci_availability_changed=all(finite(ref[f]) for f in [lower,upper])!=all(finite(rerun[f]) for f in [lower,upper]),
                ci_zero_classification_changed=excludes(ref)!=excludes(rerun) if intervals else '',
                p_availability_changed=(finite(ref.get('p_value'))!=finite(rerun.get('p_value'))) if inference_measure else False,
                bh_availability_changed=(finite(ref.get('q_primary'))!=finite(rerun.get('q_primary'))) if inference_measure else False,
                p_significance_changed=(number(ref.get('p_value'))<.05)!=(number(rerun.get('p_value'))<.05) if inference_measure and finite(ref.get('p_value')) and finite(rerun.get('p_value')) else '',
                bh_significance_changed=(number(ref.get('q_primary'))<.05)!=(number(rerun.get('q_primary'))<.05) if inference_measure and finite(ref.get('q_primary')) and finite(rerun.get('q_primary')) else '',
                valid_replicates_difference=int(rerun['valid_replicates'])-int(ref['valid_replicates']) if family=='biomarker' and estimate=='ic_years' else '',
                lower_difference=number(rerun[lower])-number(ref[lower]) if intervals else '',
                upper_difference=number(rerun[upper])-number(ref[upper]) if intervals else ''))
    return rows

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-root',type=Path,default=ROOT)
    p.add_argument('--output',type=Path,default=ROOT/'validation/output')
    p.add_argument('--reference-dir',type=Path,required=True,help='Local aggregate manuscript reference CSVs; these are not bundled')
    p.add_argument('--families',nargs='+',choices=['main','nhanes3','cause','biomarker'],default=['main','nhanes3','cause','biomarker'])
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    paths={'main':('main_cqr.csv','metal_survival_analysis/results/results_by_imputation.csv'),
      'nhanes3':('nhanes3_cqr.csv','nhanes3_analysis/results/results_by_imputation.csv'),
      'cause':('cause_cqr.csv','cause_specific/results/results_by_imputation.csv'),
      'biomarker':('biomarker.csv','biomarker_emm/output/aggregate/results.csv')}
    comparisons=[];summary=[]
    for family in args.families:
        ref,path=paths[family]
        rows=compare(family,read(args.reference_dir/ref),read(args.run_root/path))
        comparisons+=rows
        paired=[r for r in rows if r['paired_finite']]
        summary.append(dict(family=family,cells=len(rows),paired_finite=len(paired),
          availability_changes=sum(not r['same_availability'] for r in rows),
          magnitude_flags=sum(r['magnitude_flag'] is True for r in rows),sign_changes=sum(r['sign_changed'] is True for r in rows),
          ci_availability_changes=sum(r['ci_availability_changed'] for r in rows),
          ci_zero_classification_changes=sum(r['ci_zero_classification_changed'] is True for r in rows),
          p_availability_changes=sum(r['p_availability_changed'] for r in rows),
          bh_availability_changes=sum(r['bh_availability_changed'] for r in rows),
          p_significance_changes=sum(r['p_significance_changed'] is True for r in rows),
          bh_significance_changes=sum(r['bh_significance_changed'] is True for r in rows),
          replicate_count_changes=sum(r['valid_replicates_difference'] not in ['',0] for r in rows),
          maximum_absolute_difference_years=max([abs(r['difference_years']) for r in paired],default=None),
          maximum_interval_bound_difference_years=max([abs(r[k]) for r in rows
            for k in ['lower_difference','upper_difference'] if r[k]!=''],default=None)))
    output(args.output/'cell_comparison.csv',comparisons)
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
