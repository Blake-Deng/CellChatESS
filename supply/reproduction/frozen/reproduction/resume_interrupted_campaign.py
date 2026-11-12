#!/usr/bin/env python3
"""Continue the interrupted CC22 campaign without modifying its first attempt."""
from pathlib import Path
import csv, datetime, json, os, re, subprocess, sys, time
ROOT=Path(__file__).resolve().parents[1]; NEW=ROOT/'results/new'; OLD=NEW/'primary216'; RETRY=NEW/'primary216_continuation_01'; MERGED=NEW/'primary216_merged'; AUDIT=ROOT/'audit'; LOG=ROOT/'logs'; rt=json.loads((ROOT/'server_install/runtime.json').read_text())
state_path=ROOT/'server_continuation_status.json'
if state_path.exists():raise SystemExit('Continuation already exists; inspect status rather than overwrite')
state={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'prior_campaign_pid':3941749,'prior_campaign_status':'INTERRUPTED_PROCESS_MISSING','tasks':{},'note':'All first-attempt records preserved; each incomplete input is rerun as a complete pair in a fresh directory.'}
def save():
 state['updated_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();tmp=state_path.with_suffix('.tmp');tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(state_path)
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rs,fields=None):
 if not rs and fields is None:raise RuntimeError('no rows')
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rs[0]),delimiter='\t');w.writeheader();w.writerows(rs)
def run(name,cmd,required=False):
 p=LOG/('continuation_'+name+'.log')
 if p.exists():raise RuntimeError('Refusing existing log '+str(p))
 e={'status':'RUNNING','command':list(map(str,cmd)),'log':str(p.relative_to(ROOT)),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};state['tasks'][name]=e;state['current_task']=name;save()
 with p.open('w') as f:
  child=subprocess.Popen(e['command'],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT);e['pid']=child.pid;save();rc=child.wait()
 e.update(status='PASS' if rc==0 else 'FAIL',exit_code=rc,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());save()
 if rc and required:raise RuntimeError(name+' failed; see '+str(p))
 return rc==0
try:
 for p in (RETRY,MERGED):
  if p.exists():raise RuntimeError('Fresh output required: '+str(p))
 man=read(OLD/'resolved_inputs.tsv');done=[];todo=[]
 for x in man:
  b=OLD/x['experiment_id']/'block1'
  (done if all((b/f).exists() for f in ['original/measurement.tsv','ess/measurement.tsv','ess_vs_original.tsv']) else todo).append(x)
 if len(man)!=216 or len(done)+len(todo)!=216:raise RuntimeError('Unexpected input count')
 plan=AUDIT/'primary216_continuation_01_inputs.tsv'
 if plan.exists():raise RuntimeError('Refusing existing manifest')
 write(plan,todo,list(man[0]));state.update(original_complete_pairs=len(done),continuation_pairs=len(todo),original_method_results_preserved=sum(1 for p in OLD.glob('*/block1/*/measurement.tsv')));save()
 run('primary216_remaining',[sys.executable,'-u','reproduction/run_unified.py','primary216','--runtime','server_install/runtime.json','--input-manifest',plan,'--output',RETRY,'--jobs','8','--execute'],required=True)
 result=json.loads((RETRY/'completion.json').read_text())
 if result['status']!='COMPLETED' or result['method_processes_completed']!=2*len(todo):raise RuntimeError('Continuation incomplete')
 # Merge only full pairs. Partial first-attempt methods remain available in OLD.
 MERGED.mkdir(); observations=[]
 for x in man:
  source=OLD if x in done else RETRY; directory=source/x['experiment_id'];(MERGED/x['experiment_id']).symlink_to(directory,target_is_directory=True)
  for method in ['original','ess']:
   d=directory/'block1'/method; z=read(d/'measurement.tsv')[0];times=(d/'gnu_time.txt').read_text();rss=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',times);wall=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)',times)
   if not rss:raise RuntimeError('Missing RSS: '+str(d))
   z.update(experiment_id=x['experiment_id'],dataset_id=x['dataset_id'],scale=x['scale'],replicate=x['replicate'],block=1,peak_RSS_GiB=int(rss[1])/1024**2,full_process_wallclock_raw=wall[1] if wall else '',result_path=str(d/'result.rds'));observations.append(z)
 observations.sort(key=lambda z:(z['experiment_id'],z['method']));write(MERGED/'all_measurements.tsv',observations)
 for f in ['resolved_inputs.tsv','execution_plan.tsv','input_hash_verification.tsv','runtime.json']:(MERGED/f).symlink_to(OLD/f)
 original_prov=json.loads((OLD/'provenance.json').read_text());new_prov=json.loads((RETRY/'provenance.json').read_text())
 if original_prov['candidate_source_manifest_sha256']!=new_prov['candidate_source_manifest_sha256']:raise RuntimeError('Candidate source changed between attempts')
 original_prov.update(full_run=True,merge_scope='49 complete original pairs and complete continuation pairs; partial original attempts retained separately',attempts=[str(OLD),str(RETRY)],concurrency_limit=8)
 (MERGED/'provenance.json').write_text(json.dumps(original_prov,indent=2)+'\n')
 (MERGED/'completion.json').write_text(json.dumps({'status':'COMPLETED','method_processes_completed':len(observations),'full_run':True,'baseline':'CellChat 2.2.0.9001','failures':[],'attempts':[str(OLD),str(RETRY)]},indent=2)+'\n')
 state['tasks']['merge_primary216']={'status':'PASS','complete_pairs':len(man),'evidence':str((MERGED/'completion.json').relative_to(ROOT))};save()
 smoke=run('autozyme192_smoke',[sys.executable,'-u','reproduction/run_unified.py','autozyme192','--runtime','server_install/runtime.json','--output',NEW/'autozyme_smoke','--limit','1','--jobs','1','--execute'])
 if smoke:run('autozyme192_full',[sys.executable,'-u','reproduction/run_unified.py','autozyme192','--runtime','server_install/runtime.json','--output',NEW/'autozyme192','--jobs','4','--execute'])
 else:state['tasks']['autozyme192_full']={'status':'BLOCKED','reason':'AutoZyme smoke failed'};save()
 for name,script in [('ablation','run_server_ablation.py'),('xenium_target_traces','run_server_traces.py'),('xenium_full_profile','run_server_xenium.py'),('final_statistics_figures','finalize_server_results.py')]:run(name,[sys.executable,'-u',ROOT/'reproduction'/script])
 state['current_task']='FINISHED';state['status']='COMPLETE' if all(x['status']=='PASS' for x in state['tasks'].values()) else 'COMPLETE_WITH_UNRESOLVED_ITEMS';save()
except BaseException as e:
 state.update(status='STOPPED',error=repr(e));save();raise
