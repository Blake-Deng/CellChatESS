#!/usr/bin/env python3
"""New final-candidate runs. Defaults to plan-only; standard library only."""
import argparse,csv,hashlib,json,os,platform,subprocess,sys,time,threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[1]
P=argparse.ArgumentParser(description=__doc__)
P.add_argument('mode',choices=['primary216','autozyme192'])
P.add_argument('--runtime',type=Path,default=ROOT/'reproduction/config/runtime.example.json')
P.add_argument('--input-manifest',type=Path,help='Explicit manifest for a separately recorded continuation; it must retain the original hashes and rows')
P.add_argument('--path-map',type=Path)
P.add_argument('--output',type=Path,required=True)
P.add_argument('--jobs',type=int,default=1)
P.add_argument('--limit',type=int,help='Smoke only: take first N manifest inputs; does not complete full benchmark')
P.add_argument('--execute',action='store_true')
A=P.parse_args()
if A.jobs<1 or (A.limit is not None and A.limit<1):P.error('jobs/limit must be positive')
if A.jobs>(8 if A.mode=='primary216' else 4):P.error('jobs exceeds historical pair/input concurrency cap')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def readtsv(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def writetsv(p,rows):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
manifest=A.input_manifest or ROOT/'reproduction/config'/('primary_216_CC22.tsv' if A.mode=='primary216' else 'autozyme_32_CC22.tsv')
rows=readtsv(manifest)
if A.limit:rows=rows[:A.limit]
runtime=json.loads(A.runtime.read_text())
path_map=json.loads(A.path_map.read_text()) if A.path_map else {}
for r in rows:
 for old,new in sorted(path_map.items(),key=lambda x:-len(x[0])):
  if r['input_file'].startswith(old):r['input_file']=new+r['input_file'][len(old):];break
orders=[('original','autozyme','ess'),('autozyme','ess','original'),('ess','original','autozyme'),
 ('original','ess','autozyme'),('ess','autozyme','original'),('autozyme','original','ess')]
plan=[{'input_id':r['experiment_id'],'input_file':r['input_file'],'input_sha256':r['input_sha256'],
 'block':b,'order':','.join(order)} for r in rows for b,order in
 ([(1,tuple(r['method_order'].split(',')))] if A.mode=='primary216' else list(enumerate(orders,1)))]
if not A.execute:
 print(json.dumps({'status':'PLAN_ONLY','manifest_inputs':len(rows),'blocks':len(plan),
  'method_processes':len(plan)*(2 if A.mode=='primary216' else 3),'full_run':A.limit is None,'baseline':'CellChat 2.2.0.9001',
  'candidate_version':'0.2.2','first_blocks':plan[:3]},indent=2));sys.exit()
if platform.system()!='Linux':P.error('Execution requires Linux /usr/bin/time -v; planning works elsewhere')
for k in ['rscript','ess_library','cellchat_22_library']:
 if '/ABSOLUTE/' in runtime[k] or not Path(runtime[k]).exists():P.error('Set valid runtime path: '+k)
if A.mode=='autozyme192' and ('/ABSOLUTE/' in runtime['autozyme_library'] or not Path(runtime['autozyme_library']).exists()):P.error('Set autozyme_library')
if not Path('/usr/bin/time').exists():P.error('GNU time not found')
source_rows=readtsv(ROOT/'audit/candidate_source_manifest.tsv')
for r in source_rows:
 if sha(ROOT/'software/CellChatESS'/r['path'])!=r['sha256']:P.error('Candidate source changed; freeze/regenerate provenance first: '+r['path'])
source_id=sha(ROOT/'audit/candidate_source_manifest.tsv')
# Check all hashes before the first expensive experiment. Original files are never replaced.
checked={}
for r in rows:
 p=Path(r['input_file'])
 if not p.exists():P.error('Missing server input: '+str(p))
 if str(p) not in checked:checked[str(p)]=sha(p)
 if checked[str(p)]!=r['input_sha256']:P.error('Input hash mismatch: '+str(p))
out=A.output.resolve()
if out.exists():P.error('Choose a new output directory; refusing to overwrite')
out.mkdir(parents=True)
(out/'runtime.json').write_text(json.dumps(runtime,indent=2)+'\n')
writetsv(out/'resolved_inputs.tsv',rows);writetsv(out/'execution_plan.tsv',plan)
writetsv(out/'input_hash_verification.tsv',[dict(input_file=r['input_file'],expected_sha256=r['input_sha256'],actual_sha256=checked[r['input_file']],status='PASS') for r in rows])
(out/'provenance.json').write_text(json.dumps({'candidate_source_manifest_sha256':source_id,
 'candidate_version':'0.2.2','mode':A.mode,'full_run':A.limit is None and A.input_manifest is None,'baseline':'CellChat 2.2.0.9001','jobs':A.jobs,
 'hostname':platform.node(),'platform':platform.platform(),'started_unix':time.time(),
 'worker_sha256':sha(ROOT/'reproduction/inference_worker.R'),'input_manifest_sha256':sha(manifest),
 'runtime_sha256':sha(A.runtime),'driver_sha256':sha(Path(__file__)),
 'path_map':path_map},indent=2)+'\n')
for cmd,name in [(['uname','-a'],'uname.txt'),(['lscpu'],'lscpu.txt'),(['free','-b'],'memory_before.txt')]:
 try:(out/name).write_text(subprocess.run(cmd,capture_output=True,text=True,check=True).stdout)
 except (OSError,subprocess.CalledProcessError) as e:(out/name).write_text(str(e))
env=os.environ.copy()
if runtime.get('dependency_library'):env['R_LIBS']=runtime['dependency_library']+os.pathsep+env.get('R_LIBS','')
env["LC_ALL"]="C"
env['AUTOZYME_DISABLED']='true'
env['CELLCHATESS_EXPECTED_SOURCE_ID']=runtime['ess_source_id']
progress_lock=threading.Lock()
def journal(event,**kwargs):
 with progress_lock:
  with (out/'progress.jsonl').open('a') as handle:handle.write(json.dumps(dict(epoch=time.time(),event=event,**kwargs))+'\n')
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','RCPP_PARALLEL_NUM_THREADS','AUTOZYME_THREADS','ZYME_THREADS']:env[k]='1'
cc_version='2.2.0.9001'
cc_lib=runtime['cellchat_22_library']
def command(cmd,dst,name):
 (dst/(name+'.command.json')).write_text(json.dumps(cmd,ensure_ascii=False)+'\n')
 with (dst/(name+'.log')).open('w') as f:subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
def run_input(r):
 dst=out/r['experiment_id'];dst.mkdir()
 all_measurements=[]
 for b,order in ([(1,tuple(r['method_order'].split(',')))] if A.mode=='primary216' else list(enumerate(orders,1))):
  block=dst/f'block{b}';block.mkdir()
  for method in order:
   d=block/method;d.mkdir()
   cmd=['/usr/bin/time','-v','-o',str(d/'gnu_time.txt'),runtime['rscript'],'--vanilla',
    str(ROOT/'reproduction/inference_worker.R'),method,r['input_file'],str(d),cc_version,
    runtime['ess_library'],cc_lib,runtime.get('autozyme_library',cc_lib),'1',source_id,r['input_sha256']]
   journal('method_started',input_id=r['experiment_id'],block=b,method=method)
   command(cmd,d,'run')
   z=readtsv(d/'measurement.tsv')[0]
   for key in ['cells','K','nLR']:
    if int(z[key])!=int(r[key]):raise ValueError(f'{r["experiment_id"]} {key} differs from manifest')
   resource=(d/'gnu_time.txt').read_text()
   import re
   rss=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',resource)
   wall=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)',resource)
   if not rss:raise ValueError('GNU time RSS missing')
   z.update(experiment_id=r['experiment_id'],dataset_id=r['dataset_id'],scale=r['scale'],
    replicate=r['replicate'],block=b,peak_RSS_GiB=int(rss[1])*1024/1024**3,
    full_process_wallclock_raw=wall[1] if wall else '',result_path=str(d/'result.rds'))
   all_measurements.append(z)
   journal('method_completed',input_id=r['experiment_id'],block=b,method=method,inference_call_sec=z['inference_call_sec'],peak_RSS_GiB=z['peak_RSS_GiB'])
  for method in order:
   if method=='original':continue
   command([runtime['rscript'],'--vanilla',str(ROOT/'reproduction/compare_results.R'),
    str(block/'original/result.rds'),str(block/method/'result.rds'),str(block/(method+'_vs_original.tsv'))],block,method+'_comparison')
 for f in dst.rglob('*.tsv'):
  if f.name=='loaded_dlls.tsv':
   dllrows=readtsv(f)
   for z in dllrows:
    p=Path(z['path']);z['sha256']=sha(p) if p.is_file() else ''
   writetsv(f.with_name('loaded_dlls_sha256.tsv'),dllrows)
 return all_measurements
measurements=[];failures=[]
with ThreadPoolExecutor(max_workers=A.jobs) as pool:
 tasks={pool.submit(run_input,r):r for r in rows}
 for future in as_completed(tasks):
  r=tasks[future]
  try:measurements.extend(future.result());print('DONE',r['experiment_id'],flush=True)
  except Exception as e:failures.append({'input_id':r['experiment_id'],'error':repr(e)});print('FAILED',r['experiment_id'],repr(e),flush=True)
if measurements:
 measurements.sort(key=lambda x:(x['experiment_id'],int(x['block']),x['method']))
 writetsv(out/'all_measurements.tsv',measurements)
(out/'completion.json').write_text(json.dumps({'status':'FAILED' if failures else 'COMPLETED',
 'method_processes_completed':len(measurements),'full_run':A.limit is None and A.input_manifest is None,'baseline':'CellChat 2.2.0.9001','failures':failures},indent=2)+'\n')
if failures:sys.exit(1)
