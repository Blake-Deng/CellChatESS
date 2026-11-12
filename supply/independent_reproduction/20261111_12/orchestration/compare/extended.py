"""One distinct on-disk campaign for regenerated preparation or for the legacy reference."""
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import hashlib,os,sys,json,shutil,csv,argparse,subprocess,time
from statistics import median

def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def save(p,rows):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def status(path,j):
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(j,indent=2)+'\n');tmp.replace(path)
def main():
 ap=argparse.ArgumentParser(__doc__);ap.add_argument('mode',choices=['prep','legacy']);ap.add_argument('--root',type=Path,required=True);ap.add_argument('--jobs',type=int,default=1);ap.add_argument('--limit',type=int,default=0);a=ap.parse_args()
 root=a.root.resolve();here=Path(__file__).resolve().parent
 if not 1<=a.jobs<=(2 if a.mode=='prep' else 8):ap.error('concurrency cap')
 out=root/'results'/('crosscheck/preparation129' if a.mode=='prep' else 'legacy_cc16_216');out=out.with_name(out.name+('_smoke' if a.limit else ''))
 out.mkdir(exist_ok=True,parents=True)
 rr=root/'repository'/'supply';audit=rr/'reproduction/input_preparation'
 cc=rr/'legacy/CC16_ESS011'
 workers=read(root/'workspace/reproduction/config/primary_216_CC22.tsv');pmap={r['experiment_id']:r for r in workers}
 defs=read(audit/'prepared_input_audit_216.tsv')
 assert len(pmap)==216 and set(pmap)==set(r['experiment_id'] for r in defs)
 for r in defs:
  v=pmap[r['experiment_id']];r['input_file']=v['input_file']
  assert r['prepared_input_sha256']==v['input_sha256']
  assert Path(v['input_file']).is_file()
  r['source_file']=str(Path('/data/shared')/r['source_seurat_object'])
 env=os.environ.copy()
 for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','RCPP_PARALLEL_NUM_THREADS']:env[k]='1'
 env.update(LC_ALL='C',AUTOZYME_DISABLED='true')
 prefix=Path('/home/dzf/miniforge3/envs/cellchat-acceleration')
 env['PATH']=str(prefix/'bin')+os.pathsep+env['PATH']
 libbase=prefix/'lib/R/library'
 # Preserve the actual July CellChat installation for historical preparation and old timings.
 env['R_LIBS_USER']=str(libbase);env['R_LIBS']=str(libbase)
 rs=str(prefix/'bin/Rscript');cmd=[rs,'--vanilla']
 basecheck="stopifnot(as.character(packageVersion('CellChat'))=='1.6.1',as.character(packageVersion('Seurat'))=='5.5.1',normalizePath(find.package('CellChat'))==normalizePath(file.path(Sys.getenv('R_LIBS'),'CellChat'))); print(sessionInfo()); print(packageDescription('CellChat')[c('Version','RemoteSha','Built')]);"
 with (out/'preflight.txt').open('w') as f:subprocess.run(cmd+['-e',basecheck],env=env,stdout=f,stderr=f,check=True)
 archive=here/'CellChatESS_0.1.1.tar.gz'
 minimum_disk=(22 if a.mode=='prep' else 10)*2**30
 available_kb=int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:')))
 if shutil.disk_usage(root).free<minimum_disk or available_kb/1024**2<(24 if a.mode=='prep' else 6)*a.jobs:
  raise RuntimeError(f'Insufficient free disk/RAM for {a.mode} with {a.jobs} jobs')
 if a.mode=='legacy':
  assert sha(archive)=='603b81faf93abc3b5ef6b85b169a0f96bd659da1915be25e1cd18e080f3e1f48'
  lib=root/'libraries/legacy_cc16_ess011';lib.mkdir(parents=True,exist_ok=True)
  with (out/'install_ess011.log').open('w') as f:subprocess.run([str(prefix/'bin/R'),'CMD','INSTALL','--preclean','--clean','--library='+str(lib),str(archive)],env=env,stdout=f,stderr=f,check=True)
  env['CELLCHAT_MAIN_R_LIB']=str(lib)
  chosen=sorted(defs,key=lambda r:r['experiment_id'])
  if a.limit:chosen=chosen[:a.limit]
 if a.mode=='prep':
  ranks=lambda r:(r['target_size']!='all',int(r['repeat']),r['experiment_id'])
  canon={}
  for r in sorted(defs,key=ranks):canon.setdefault(r['prepared_input_sha256'],r)
  assert len(canon)==129
  save(out/'all216_to_129.tsv',[dict(r,regeneration_id=canon[r['prepared_input_sha256']]['experiment_id'],regeneration_seed=canon[r['prepared_input_sha256']]['sampling_seed']) for r in defs])
  chosen=sorted(canon.values(),key=lambda r:(int(r['actual_cells']),r['experiment_id']))
  if a.limit:chosen=chosen[:a.limit]
  sources=sorted(set(Path(r['source_file']) for r in chosen));save(out/'source_inputs.tsv',[dict(path=str(p),size_bytes=p.stat().st_size,current_sha256=sha(p)) for p in sources])
 save(out/'plan.tsv',chosen)
 count=0;failures=[];rows=[];start=time.time()
 def process(r):
  eid=r['experiment_id'];target=out/'records'/eid;target.mkdir(parents=True,exist_ok=True)
  result=target/'verified.json'
  if result.is_file():return json.loads(result.read_text())
  if any(target.iterdir()):
   target.rename(target.with_name(target.name+'_prior_attempt_'+str(time.time_ns())))
   target.mkdir()
  frozen=Path(r['input_file']);assert sha(frozen)==r['prepared_input_sha256']
  def exec_(args,name):
   with (target/(name+'.log')).open('w') as log:subprocess.run(args,env=env,stdout=log,stderr=log,check=True)
  if a.mode=='prep':
   exec_(cmd+[str(here/'prep_one.R'),str(audit/'cellchat_benchmark.R'),r['source_file'],r['target_size'],r['sampling_seed'],r['input_file'],str(target),eid],'prepare')
   tests=read(target/'prep_validation.tsv');same=all((v['shape_names'] in ('1','TRUE') and float(v['n_diff'])==0 and float(v['max_abs'])==0) if v['slot'] in ('data','data.signaling') else v['identical'] in ('1','TRUE') for v in tests)
   data=dict(id=eid,status='PASS' if same else 'DIFFERENT',branch=tests[0]['branch'],cells=int(r['actual_cells']),expected_sha256=r['prepared_input_sha256'],fresh_serialized_sha256=sha(target/'regenerated_cellchat.rds'),group_and_array_slot_agreement=same)
  else:
   ref=frozen.with_name('baseline.computed.rds');assert ref.is_file()
   obs={}
   for method in ('original','ess'):
    t=target/(method+'.tsv')
    exec_(['/usr/bin/time','-v','-o',str(target/(method+'.time'))]+cmd+[str(cc/'scripts/common_timer_worker.R'),'--method',method,'--id',eid,'--input',str(frozen),'--reference',str(ref),'--out',str(t),'--result',str(target/(method+'.rds'))],method)
    rr_=read(t)[0];assert rr_['CellChat_version']=='1.6.1' and rr_['CellChatESS_version']=='0.1.1';obs[method]=float(rr_['inference_call_sec'])
   exec_(cmd+[str(here/'legacy_count.R'),str(target/'original.rds'),str(target/'ess.rds'),str(target/'outputs.tsv')],'compare')
   old=Path('/data/dzf/cellchatess_revision_20261026/results/common_external_timer_20261029/full/pairs')/eid
   previous_ok=True;previous_checked=0
   for method in ('original','ess'):
    archived=old/(method+'.rds')
    if not archived.is_file():continue
    exec_(cmd+[str(here/'legacy_count.R'),str(archived),str(target/(method+'.rds')),str(target/(method+'_vs_archived.tsv'))],method+'_archive_check')
    direct=read(target/(method+'_vs_archived.tsv'))
    previous_ok=previous_ok and all(z['dimensions_names_finite']=='TRUE' and int(z['changed'])==0 for z in direct)
    previous_checked+=1
   tests=read(target/'outputs.tsv');prob=float(tests[0]['max_abs']);same=all(v['dimensions_names_finite']=='TRUE' for v in tests) and all(int(v['changed'])==0 for v in tests if v['field']!='prob') and prob<=1.395e-16 and previous_ok
   data=dict(id=eid,dataset=r['dataset_id'],status='PASS' if same else 'DIFFERENT',original_sec=obs['original'],ess_sec=obs['ess'],speedup=obs['original']/obs['ess'],max_probability_difference=prob,archived_method_files_compared=previous_checked,all_available_archived_values_identical=previous_ok,changed_pvalue_elements=tests[1]['changed'],changed_count_elements=tests[2]['changed'])
  status(result,data);return data
 with ThreadPoolExecutor(max_workers=a.jobs) as ex:
  fs={ex.submit(process,r):r for r in chosen}
  for future in as_completed(fs):
   rid=fs[future]['experiment_id']
   try:rows.append(future.result())
   except Exception as e:failures.append(dict(id=rid,error=repr(e)))
   count+=1;status(out/'progress.json',dict(done=count,planned=len(chosen),exceptions=failures,observed_mismatches=sum(v['status']!='PASS' for v in rows),elapsed_seconds=time.time()-start))
 if rows:save(out/'results.tsv',sorted(rows,key=lambda r:r['id']))
 passed=len(rows)==len(chosen) and all(v['status']=='PASS' for v in rows)
 report=dict(status='PASS' if passed else 'FAIL',planned=len(chosen),processed=len(rows),failures=failures)
 if rows and a.mode=='legacy':
  import numpy as np
  groups=sorted({r['dataset'] for r in rows});g=[np.array([r['speedup'] for r in rows if r['dataset']==dataset]) for dataset in groups];rng=np.random.default_rng(20261026)
  med=[np.median(np.concatenate([g[i] for i in rng.integers(0,len(g),len(g))])) for rep in range(10000)]
  report.update(pooled_median=median(r['speedup'] for r in rows),dataset_cluster_percentile_ci=np.percentile(med,[2.5,97.5]).tolist(),max_probability_difference=max(r['max_probability_difference'] for r in rows),archived_median=70.82,reference_version='Historical installed CellChat 1.6.1, July build, upstream Git commit not recorded; archived ESS 0.1.1 release tarball, compiled in new isolated library.')
 status(out/'completion.json',report);return int(not passed)
if __name__=='__main__':sys.exit(main())
