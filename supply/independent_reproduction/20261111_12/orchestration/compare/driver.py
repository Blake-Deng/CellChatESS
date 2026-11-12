"""Run independent validation phases in sequence even if an earlier phase reports a discrepancy."""
import argparse,json,sys,subprocess,os,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);a=p.parse_args();r=a.root.resolve();src=Path(__file__).resolve().parent;new=r/'workspace/results/new';o=r/'results/crosscheck';o.mkdir(exist_ok=True,parents=True);rt=json.loads((r/'workspace/server_install/runtime.json').read_text());py='/home/dzf/miniforge3/bin/python';basic='/usr/bin/python3'
common=['--root',str(r)];cmds=[
 ('record_and_array_comparison',[basic,str(src/'core_audit.py'),'--repo',str(r/'repository'),'--primary',str(new/'primary216'),'--autozyme',str(new/'autozyme192'),'--out',str(o/'core'),'--rscript',rt['rscript']]),
 ('fresh_tables_figures_and_uncertainty',[py,str(src/'fresh_figures.py'),*common]),
 ('historical_CC16_ESS011_216',[basic,str(src/'extended.py'),'legacy',*common,'--jobs','8']),
 ('regenerate_129_prepared_inputs',[basic,str(src/'extended.py'),'prep',*common,'--jobs','2'])]
env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','BLIS_NUM_THREADS':'1','MPLCONFIGDIR':str(src/'mpl_cache'),'PYTHONUNBUFFERED':'1'}
res=[];trace=o/'phase_log.json'
def save():
 dst=trace.with_suffix('.json.tmp');dst.write_text(json.dumps(res,indent=2)+'\n');dst.replace(trace)
last=json.loads(trace.read_text()) if trace.exists() else []
for name,cmd in cmds:
 # A restart reuses only explicitly completed phases whose products survive.
 manifests={'record_and_array_comparison':o/'core/audit.json','fresh_tables_figures_and_uncertainty':r/'results/figures/FigureS1_fresh.pdf','historical_CC16_ESS011_216':r/'results/legacy_cc16_216/completion.json','regenerate_129_prepared_inputs':o/'preparation129/completion.json'}
 prev=next((x for x in last if x['stage']==name and x.get('exit_code')==0),None)
 if prev and manifests[name].is_file():res.append(prev);continue
 reject_dir={'record_and_array_comparison':o/'core','fresh_tables_figures_and_uncertainty':r/'results/figures'}.get(name)
 if reject_dir is not None and reject_dir.exists():reject_dir.rename(reject_dir.with_name(reject_dir.name+'_prior_attempt_'+str(time.time_ns())))
 dest=o/(name+'_'+str(time.time_ns())+'.log')
 started=time.time();rrecord={'stage':name,'argv':cmd,'started':started,'log':str(dest),'exit_code':None};res.append(rrecord)
 (o/'current_stage.json').write_text(json.dumps(rrecord,indent=2))
 with dest.open('w') as f:pc=subprocess.run(cmd,env=env,stdout=f,stderr=f)
 rrecord.update(exit_code=pc.returncode,elapsed_s=time.time()-started);save();print(name,pc.returncode,flush=True)
save()
final=subprocess.run([basic,str(src/'report.py'),*common],env=env)
sys.exit(final.returncode)
