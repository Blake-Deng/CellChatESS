"""Audit fresh core blocks against the pinned reference. Read-only on both inputs."""
import argparse,csv,hashlib,json,math,subprocess,sys
from collections import Counter
from pathlib import Path
from statistics import median

def tab(path):
    with path.open() as f:
        return list(csv.DictReader(f,delimiter='\t'))
def save(path,rows):
    if not rows:return
    with path.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def index(rows,keys):
    ix={tuple(r[k] for k in keys):r for r in rows}
    if len(ix)!=len(rows):raise ValueError(f'duplicate keys {keys}')
    return ix
def truth(s):return str(s).lower()=='true'
def equal(x,y):
    if x in ('','NA','NaN') and y in ('','NA','NaN'):return True
    if x==y:return True
    try:return float(x)==float(y)
    except (ValueError,TypeError):return False
def read_json(p):return json.loads(p.read_text())
def run():
    ap=argparse.ArgumentParser(__doc__);ap.add_argument('--repo',required=True,type=Path);ap.add_argument('--primary',type=Path,required=True);ap.add_argument('--autozyme',type=Path,required=True);ap.add_argument('--out',required=True,type=Path);ap.add_argument('--rscript',required=True);a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=False)
    s=a.repo/'supply'; refs=tab(s/'data/all_main_measurements_1008.tsv'); compar=tab(s/'data/numerical_comparisons.tsv')
    perf=[];rdspairs=[];archdrift=[];missing=[];fail=[];reports={};fold={}
    for queue,n,root,want in [('primary216','primary',a.primary,432),('autozyme192','autozyme',a.autozyme,576)]:
        base=s/'results'/('primary216_merged' if queue=='primary216' else 'autozyme192')
        complete=read_json(root/'completion.json')
        if complete.get('status')!='COMPLETED' or not complete.get('full_run') or complete['method_processes_completed']!=want:raise ValueError(f'incomplete {queue}')
        ref=[r for r in refs if r['queue']==queue]
        if len(ref)!=want:raise ValueError(f'reference count {queue}')
        k=['experiment_id','block','method']
        x=index(ref,k);y=index(tab(root/'all_measurements.tsv'),k)
        if x.keys()!=y.keys():raise ValueError(f'wrong measurement set {queue}')
        order=lambda p:{(r['input_id'],r['block'],r['order'],r['input_sha256']) for r in tab(p/'execution_plan.tsv')}
        if order(root)!=order(base):raise ValueError(f'method order or hashes changed {queue}')
        hp=tab(root/'input_hash_verification.tsv')
        if len(hp)!=(216 if n=='primary' else 32) or not all(r['expected_sha256']==r['actual_sha256'] and r['status']=='PASS' for r in hp):raise ValueError('hash check failed')
        expected=set(r['input_sha256'] for r in ref)
        if set(r['actual_sha256'] for r in hp)!=expected:raise ValueError('hashed inputs differ')
        for key,r in x.items():
            v=y[key]
            for field in ['input_sha256','candidate_source_id','CellChat_version','ESS_loaded_version','AutoZyme_loaded_version','cells','signaling_genes','K','nLR','probability_payload_bytes','pvalue_payload_bytes','nboot','seed_use','threads_requested','threads_used','timing_scope','dataset_id','scale','replicate']:
                if not equal(r[field],v[field]):fail.append(dict(queue=queue,record=str(key),field=field,reference=r[field],fresh=v[field]))
            eid,b,m=key
            old=Path(r['result_path']);new=root/eid/('block'+b)/m/'result.rds'
            if not new.is_file():raise FileNotFoundError(new)
            pair=dict(queue=queue,experiment_id=eid,block=b,method=m,archive_file=str(old),fresh_file=str(new))
            rdspairs.append(pair)
            if not old.is_file():missing.append(pair)
            t,t0,mem,mem0=[float(w) for w in (v['inference_call_sec'],r['inference_call_sec'],v['peak_RSS_GiB'],r['peak_RSS_GiB'])]
            if not all(math.isfinite(val) and val>0 for val in [t,t0,mem,mem0]):raise ValueError('invalid time/RSS')
            perf.append(dict(queue=queue,experiment_id=eid,block=b,method=m,new_inference_s=t,archive_inference_s=t0,new_to_archive_time=t/t0,new_peak_rss_GiB=mem,archive_peak_rss_GiB=mem0,new_to_archive_rss=mem/mem0))
            fold.setdefault((queue,eid,b),{})[m]=(t,t0,mem,mem0)
        fresh=[]
        for path in sorted(root.glob('*/block*/*_vs_original.tsv')):
            for v in tab(path):fresh.append(dict(v,input_id=path.parent.parent.name,block=path.parent.name,comparison=path.stem))
        x=index([r for r in compar if r['queue']==n],['input_id','block','comparison','output'])
        y=index(fresh,['input_id','block','comparison','output'])
        if x.keys()!=y.keys():raise ValueError(f'wrong per-block comparator set {queue}')
        pe=[]
        for key,r in x.items():
            v=y[key];field=key[-1]
            if not all(truth(v[k]) for k in ['dimensions_equal','dimnames_equal','finite_values']):fail.append(dict(queue=queue,record=str(key),field='dimensions_names_finite',reference='all true',fresh=str(v)))
            if field in ('pval','count') and (int(v['different_elements'])!=0 or not truth(v['elementwise_exact'])):fail.append(dict(queue=queue,record=str(key),field=field,reference='exact',fresh=v['different_elements']))
            if field=='prob':pe.append(float(v['max_absolute_difference']))
            for f in ['max_absolute_difference','different_elements','elementwise_exact']:
                same=(truth(v[f])==truth(r[f])) if f=='elementwise_exact' else equal(v[f],r[f])
                if not same:archdrift.append(dict(queue=queue,record=str(key),field=f,reference=r[f],fresh=v[f]))
        reports[queue]=dict(completed_methods=want,paired_output_checks=len(y),max_probability_error=max(pe),distinct_input_sha256=len(expected))
        if max(pe)>1.11e-16+0.005e-16:fail.append(dict(queue=queue,record='pooled',field='max_probability_difference',reference='paper rounded 1.11e-16',fresh=str(max(pe))))
    save(a.out/'paired_method_files.tsv',rdspairs)
    save(a.out/'unavailable_archived_method_files.tsv',missing)
    save(a.out/'timing_rss_observations.tsv',perf)
    save(a.out/'fresh_within_run_failures.tsv',fail)
    save(a.out/'change_in_archived_summaries.tsv',archdrift)
    direct=a.out/'direct_array_checks.tsv'
    if len(rdspairs)>len(missing):subprocess.run([a.rscript,'--vanilla',str(Path(__file__).with_name('arrays.R')),str(a.out/'paired_method_files.tsv'),str(direct)],check=True)
    audit=tab(direct) if direct.exists() else []
    if len(audit)!=(len(rdspairs)-len(missing))*3:raise ValueError('wrong direct comparison count')
    wrong=[r for r in audit if not truth(r['value_identical'])]
    save(a.out/'changed_archived_arrays.tsv',wrong)
    for q in reports:
        name='autozyme' if q=='autozyme192' else 'original'
        vals=[d for (qq,e,b),d in fold.items() if qq==q]
        reports[q].update(fresh_median_comparator_over_ess=median(d[name][0]/d['ess'][0] for d in vals),archive_median_comparator_over_ess=median(d[name][1]/d['ess'][1] for d in vals),fresh_ess_faster=sum(d[name][0]>d['ess'][0] for d in vals),fresh_ess_lower_RSS=sum(d[name][2]>d['ess'][2] for d in vals))
    j={'status':'MISMATCH' if fail or wrong or archdrift else ('PASS_WITHOUT_ALL_ARCHIVED_ARRAYS' if missing else 'PASS'),'queues':reports,'archived_arrays_checked':len(audit),'missing_archived_method_files':len(missing),'different_archived_arrays':len(wrong),'changed_recorded_summary_rows':len(archdrift),'fresh_acceptance_failures':len(fail),'performance_comparison':'Measured again; raw time and RSS are observations and not exact-equality acceptance gates.'}
    (a.out/'audit.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j,indent=2))
    return int(bool(fail or wrong or archdrift))
if __name__=='__main__':sys.exit(run())
