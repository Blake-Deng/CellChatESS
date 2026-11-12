#!/usr/bin/env python3
"""Sequential, independent-process Original-vs-final-ESS Xenium run with safety guards."""
from __future__ import annotations
import argparse, csv, hashlib, json, os, signal, subprocess, sys, time, shutil, tarfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RUNTIME=json.loads((ROOT/'server_install/runtime.json').read_text())
INPUT=Path('/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/prepared_xenium_ovary_full_grid50.rds')
EXPECTED_INPUT_SHA='cc6155a5de127ba34b8a26e7f69e7bffc704448f976107253711623d2a88de37'
PBMC=ROOT/'validation/local_smoke/data/pbmc3k_prepared_cellchat.rds'
PBMC_SHA='1305076852f518a6cb57064e8f3678ea3b6c902290206a2463e1d07f2ca817c2'
ESS_DLL_EXPECTED='c46577016d28829e5dcd1c694b3f185a762a6f925bf05cfe0721981030085752'
MIN_START_GIB=350.0
ABORT_AVAILABLE_GIB=180.0
STOP_SIGNAL=None
def request_stop(signum, frame):
    global STOP_SIGNAL
    STOP_SIGNAL=signum

def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def available_gib()->float:
    line=next(x for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:'))
    return int(line.split()[1])/1024**2
def descendants(pid:int)->list[int]:
    pending=[pid];found=[]
    while pending:
        cur=pending.pop();task=Path('/proc')/str(cur)/'task'/str(cur)/'children'
        try: children=[int(x) for x in task.read_text().split()]
        except OSError: continue
        pending.extend(children)
        try:
            cmd=(Path('/proc')/str(cur)/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
            try:comm=(Path('/proc')/str(cur)/'comm').read_text().strip()
            except OSError:comm=''
            if 'run_xenium_paired_worker.R' in cmd and comm in ('R','Rscript'):found.append(cur)
        except OSError:pass
    return found
def snapshot(out:Path,tag:str)->None:
    d=out/'environment';d.mkdir(parents=True,exist_ok=True)
    for name,cmd in [('lscpu.txt',['lscpu']),('uname.txt',['uname','-a']),('free.txt',['free','-h']),
                     ('df_data.txt',['df','-h','/data']),('df_home.txt',['df','-h','/home']),
                     ('uptime.txt',['uptime'])]:
        with (d/f'{tag}_{name}').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=False)
    rows=[]
    for line in subprocess.run(['ps','-eo','pid,rss,comm','--sort=-rss'],capture_output=True,text=True).stdout.splitlines()[1:]:
        z=line.strip().split(None,2)
        if len(z)==3 and int(z[1])>=10*1024*1024:rows.append({'pid':z[0],'rss_gib':f'{int(z[1])/1024**2:.3f}','comm':z[2]})
    with (d/f'{tag}_other_processes_over_10GiB.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['pid','rss_gib','comm'],delimiter='\t');w.writeheader();w.writerows(rows)
    (d/f'{tag}_resources.json').write_text(json.dumps({'timestamp_epoch':time.time(),'MemAvailable_GiB':available_gib(),
      'data_free_bytes':__import__('shutil').disk_usage('/data').free,'home_free_bytes':__import__('shutil').disk_usage('/home').free,
      'unrelated_processes_over_10GiB':rows,'policy':f'start only at >= {MIN_START_GIB} GiB available; stop this job below {ABORT_AVAILABLE_GIB} GiB'},indent=2)+'\n')
def env_for()->dict[str,str]:
    e=os.environ.copy()
    e.update({'R_LIBS':RUNTIME['dependency_library'],'R_LIBS_USER':'','R_LIBS_SITE':'',
      'R_PROFILE_USER':'/dev/null','R_ENVIRON_USER':'/dev/null','AUTOZYME_DISABLED':'true',
      'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','BLIS_NUM_THREADS':'1',
      'VECLIB_MAXIMUM_THREADS':'1','RCPP_PARALLEL_NUM_THREADS':'1','AUTOZYME_THREADS':'1','ZYME_THREADS':'1',
      'LC_ALL':'C','TZ':'Asia/Shanghai','EXPECTED_ESS_DLL_SHA256':ESS_DLL_EXPECTED})
    return e
def run_worker(method:str,input_path:Path,input_sha:str,out:Path,full:bool,abort_floor:float)->dict:
    if out.exists():raise RuntimeError(f'refusing to overwrite {out}')
    out.mkdir(parents=True)
    worker=ROOT/'reproduction/run_xenium_paired_worker.R'
    cc_source=ROOT/'server_install/extracted/CellChat_2.2.0.9001/CellChat-75253cd0c9e68410e6e721a6d3a0419a1d7e358f/R/modeling.R'
    ess_source=ROOT/'software/CellChatESS/R/accel_cellchat.R'
    cmd=[RUNTIME['rscript'],'--vanilla',str(worker),method,str(input_path),input_sha,str(out),
      RUNTIME['cellchat_22_library'],RUNTIME['ess_library'],RUNTIME['dependency_library'],str(cc_source),str(ess_source)]
    env=env_for();env['REQUIRE_XENIUM_FULL']='1' if full else '0'
    wrapped=['/usr/bin/time','-v','-o',str(out/'gnu_time.txt')]+cmd
    (out/'command.json').write_text(json.dumps(wrapped,indent=2)+'\n')
    (out/'launch_environment.json').write_text(json.dumps({'method':method,'cwd':str(ROOT),'input_sha256_expected':input_sha,
      'MemAvailable_GiB_at_start':available_gib(),'safety_stop_floor_GiB':abort_floor,
      'threads':{'OMP_NUM_THREADS':1,'OPENBLAS_NUM_THREADS':1,'MKL_NUM_THREADS':1,'BLIS_NUM_THREADS':1,
       'VECLIB_MAXIMUM_THREADS':1,'RCPP_PARALLEL_NUM_THREADS':1,'ESS_threads':1,'future_workers':1},
      'AutoZyme':'not on R library path; AUTOZYME_DISABLED=true'},indent=2)+'\n')
    log=(out/'run.log').open('w')
    start=time.monotonic();p=subprocess.Popen(wrapped,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    rpid=None;sampler=None;sampler_log=None;stop_reason=None;last=time.monotonic()
    while p.poll() is None:
        if rpid is None:
            found=descendants(p.pid)
            if found:
                rpid=found[0]
                sampler_log=(out/'sampler.log').open('w')
                sampler=subprocess.Popen([sys.executable,str(ROOT/'reproduction/sample_process_resources.py'),
                  '--pid',str(rpid),'--output',str(out/'rss_samples.tsv'),'--stages',str(out/'stage_events.tsv'),'--interval','1'],
                  cwd=ROOT,env=env,stdout=sampler_log,stderr=subprocess.STDOUT,start_new_session=True)
                (out/'actual_process.json').write_text(json.dumps({'gnu_time_wrapper_pid':p.pid,'R_pid':rpid},indent=2)+'\n')
        avail=available_gib()
        if avail<abort_floor or STOP_SIGNAL is not None:
            stop_reason=(f'Controller received signal {STOP_SIGNAL}' if STOP_SIGNAL is not None else
              f'MemAvailable fell below safety floor {abort_floor:.1f} GiB (observed {avail:.1f} GiB)')
            try:os.killpg(p.pid,signal.SIGTERM)
            except ProcessLookupError:pass
            try:p.wait(timeout=60)
            except subprocess.TimeoutExpired:
                try:os.killpg(p.pid,signal.SIGKILL)
                except ProcessLookupError:pass
                p.wait()
            break
        if time.monotonic()-last>=60:
            print(f'{method}: elapsed={time.monotonic()-start:.0f}s, MemAvailable={avail:.1f} GiB, R_pid={rpid}',flush=True)
            last=time.monotonic()
        time.sleep(5)
    rc=p.wait();log.close()
    if sampler is not None:
        try:sampler.wait(timeout=70)
        except subprocess.TimeoutExpired:sampler.terminate();sampler.wait()
        if sampler_log:sampler_log.close()
    state={'method':method,'status':'COMPLETED' if rc==0 and not stop_reason else ('INTERRUPTED' if STOP_SIGNAL is not None else ('ABORTED_RESOURCE_GUARD' if stop_reason else 'FAILED')),
      'exit_code':rc,'safety_stop_reason':stop_reason,'elapsed_sec':time.monotonic()-start,'R_pid':rpid,
      'MemAvailable_GiB_at_end':available_gib(),'result_rds_exists':(out/'result.rds').is_file()}
    (out/'run_status.json').write_text(json.dumps(state,indent=2)+'\n')
    return state
def verify_static(out:Path)->dict:
    src=json.loads((ROOT/'audit/server_source_identity.json').read_text())
    ess_dll=ROOT/'server_install/library/ess022/CellChatESS/libs/CellChatESS.so'
    ccsource=ROOT/'server_install/extracted/CellChat_2.2.0.9001/CellChat-75253cd0c9e68410e6e721a6d3a0419a1d7e358f/R/modeling.R'
    repo=ROOT/'software/CellChatESS_candidate_git'
    checks={'CellChat_commit':src['upstream_commits']['CellChat_2.2.0.9001_provenance']['commit'],
      'CellChat_archive_sha256':src['upstream_commits']['CellChat_2.2.0.9001_provenance']['sha256'],
      'CellChat_modeling_R_sha256':sha(ccsource),'ESS_commit':src['git_commit'],'ESS_tree':src['git_tree'],
      'ESS_source_manifest_sha256':src['candidate_source_manifest_sha256'],'ESS_tarball_sha256':src['source_tarball_sha256'],
      'ESS_installed_DLL_path':str(ess_dll),'ESS_installed_DLL_sha256':sha(ess_dll),
      'ESS_expected_DLL_sha256':ESS_DLL_EXPECTED,'ESS_git_HEAD':subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip(),
      'ESS_git_worktree_status':subprocess.check_output(['git','-C',str(repo),'status','--porcelain'],text=True).strip(),
      'CellChat_installed_library':RUNTIME['cellchat_22_library'],'ESS_installed_library':RUNTIME['ess_library']}
    checks['ESS_installed_DLL_matches_expected']=checks['ESS_installed_DLL_sha256']==ESS_DLL_EXPECTED
    checks['ESS_git_HEAD_matches_audit']=checks['ESS_git_HEAD']==checks['ESS_commit']
    checks['ESS_git_worktree_clean']=not bool(checks['ESS_git_worktree_status'])
    manifest=ROOT/'audit/candidate_source_manifest.tsv'
    checks['ESS_source_manifest_recomputed_sha256']=sha(manifest)
    checks['ESS_source_files_revalidated']=[]
    with manifest.open(newline='') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            p=ROOT/'software/CellChatESS'/row['path']
            ok=p.is_file() and sha(p)==row['sha256'] and p.stat().st_size==int(row['size_bytes'])
            checks['ESS_source_files_revalidated'].append({'file':row['path'],'matches':ok})
    tar=ROOT/'software/CellChatESS_0.2.2.tar.gz'
    checks['ESS_tarball_recomputed_sha256']=sha(tar)
    archive=ROOT/'upstream/CellChat-75253cd0c9e68410e6e721a6d3a0419a1d7e358f.tar.gz'
    checks['CellChat_archive_recomputed_sha256']=sha(archive)
    with tarfile.open(archive,'r:gz') as tf:
        mismatches=[];checked=0
        for member in tf.getmembers():
            rel=Path(*Path(member.name).parts[1:])
            if member.isfile() and len(rel.parts)>1 and rel.parts[0] in ('R','src'):
                p=ccsource.parents[1]/rel
                if not p.is_file() or hashlib.sha256(tf.extractfile(member).read()).hexdigest()!=sha(p):mismatches.append(str(rel))
                checked+=1
    checks['CellChat_extracted_R_src_files_checked']=checked
    checks['CellChat_extracted_R_src_archive_mismatches']=mismatches
    checks['all_revalidated_source_and_archive_hashes_match']=(
      checks['ESS_source_manifest_recomputed_sha256']==checks['ESS_source_manifest_sha256'] and
      checks['ESS_tarball_recomputed_sha256']==checks['ESS_tarball_sha256'] and
      checks['CellChat_archive_recomputed_sha256']==checks['CellChat_archive_sha256'] and
      all(x['matches'] for x in checks['ESS_source_files_revalidated']) and not mismatches and checked>0)
    (out/'source_and_install_identity.json').write_text(json.dumps(checks,indent=2)+'\n')
    if not (checks['ESS_installed_DLL_matches_expected'] and checks['ESS_git_HEAD_matches_audit'] and checks['ESS_git_worktree_clean'] and checks['all_revalidated_source_and_archive_hashes_match']):
        raise RuntimeError('Static source/install identity verification failed; see source_and_install_identity.json')
    return checks
def run_compare(out:Path,expected_lr:int)->int:
    compare=ROOT/'reproduction/compare_xenium_paired_full627.R'
    targets=ROOT/'reproduction/config/xenium_five_targets.tsv'
    cmd=[RUNTIME['rscript'],'--vanilla',str(compare),str(out/'original/result.rds'),str(out/'ess/result.rds'),
      str(targets),str(out/'comparison'),str(expected_lr)]
    (out/'comparison').mkdir()
    (out/'comparison_command.json').write_text(json.dumps(cmd,indent=2)+'\n')
    env=env_for()
    with (out/'comparison.log').open('w') as f:rc=subprocess.run(cmd,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
    (out/'comparison_status.json').write_text(json.dumps({'status':'COMPLETED' if rc==0 else 'FAILED','exit_code':rc},indent=2)+'\n')
    return rc
def main()->int:
    signal.signal(signal.SIGTERM,request_stop);signal.signal(signal.SIGINT,request_stop)
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,default=ROOT/'results/new/xenium_paired_full627_20261101_strict01')
    ap.add_argument('--smoke-only',action='store_true',help='run only the paired small-input smoke')
    ap.add_argument('--resume-after-smoke',action='store_true',help='reuse completed paired smoke outputs in this directory and start full runs')
    ap.add_argument('--abort-floor-gib',type=float,default=ABORT_AVAILABLE_GIB)
    a=ap.parse_args();out=a.output.resolve()
    if out.exists() and not a.resume_after_smoke:raise SystemExit(f'Fresh output directory required; refusing to overwrite: {out}')
    if a.resume_after_smoke and not out.is_dir():raise SystemExit(f'--resume-after-smoke requires an existing output directory: {out}')
    if not a.smoke_only and available_gib()<MIN_START_GIB:
        raise SystemExit(f'BLOCKED before full run: MemAvailable={available_gib():.1f} GiB < {MIN_START_GIB:.1f} GiB')
    if not a.resume_after_smoke:out.mkdir(parents=True)
    if shutil.disk_usage(out).free<200*1024**3:raise SystemExit('BLOCKED: less than 200 GiB free disk for results')
    snapshot(out,'before_full' if a.resume_after_smoke else 'before')
    verify_static(out)
    if sha(PBMC)!=PBMC_SHA:raise SystemExit('PBMC smoke input hash mismatch')
    if not a.smoke_only and sha(INPUT)!=EXPECTED_INPUT_SHA:raise SystemExit('Full Xenium input hash mismatch')
    manifest={'created_epoch':time.time(),'scope':'strict paired Original CellChat 2.2.0.9001 vs final CellChatESS 0.2.2',
      'input_path':str(INPUT),'expected_input_sha256':EXPECTED_INPUT_SHA,'smoke_input_path':str(PBMC),'smoke_input_sha256':PBMC_SHA,
      'parameters':{'type':'triMean','raw.use':True,'population.size':False,'nboot':100,'seed.use':1,'Kh':0.5,'n':1,
        'ESS_algorithm':'ess','ESS_threads':1,'future_workers':1},
      'method_order':['original','ess'],'independent_R_processes':True,'concurrent_method_processes':1,
      'same_R_read_save_worker':str(ROOT/'reproduction/run_xenium_paired_worker.R'),
      'serialization':{'compress':'gzip','version':3,'fields':['prob','pval','count','groups','lr_names','parameters']},
      'timer_definitions':{'inference':'proc.time elapsed around requested computeCommunProb call; excludes startup, package load, input read, aggregateNet and serialization',
        'full_process':'GNU time -v around Rscript from process launch through input read, inference, aggregateNet, serialization and process exit'},
      'input_sha256_verified_before_and_inside_each_R_process':True,'AutoZyme':'not activated and its R library is excluded'}
    (out/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if available_gib()<MIN_START_GIB and not a.smoke_only:
        raise SystemExit(f'BLOCKED before smoke/full run after environment snapshot: MemAvailable={available_gib():.1f} GiB')
    # Small paired smoke runs exercise the exact worker, function dispatch, save format and comparator.
    if not a.resume_after_smoke:
        for method in ('original','ess'):
            d=out/'smoke'/method
            if not d.parent.exists():d.parent.mkdir(parents=True)
            state=run_worker(method,PBMC,PBMC_SHA,d,False,a.abort_floor_gib)
            if state['status']!='COMPLETED':
                (out/'smoke_status.json').write_text(json.dumps({'status':'FAILED','failed_method':method,**state},indent=2)+'\n')
                return 2
    else:
        for method in ('original','ess'):
            state_path=out/'smoke'/method/'run_status.json'
            result_path=out/'smoke'/method/'result.rds'
            if not state_path.is_file() or not result_path.is_file() or json.loads(state_path.read_text()).get('status')!='COMPLETED':
                raise SystemExit('Cannot resume: completed smoke output missing or failed for '+method)
    smoke_compare=ROOT/'reproduction/compare_xenium_paired_full627.R'
    smoke_out=out/'smoke/comparison_retry01' if a.resume_after_smoke else out/'smoke/comparison'
    smoke_out.mkdir()
    smoke_cmd=[RUNTIME['rscript'],'--vanilla',str(smoke_compare),str(out/'smoke/original/result.rds'),
      str(out/'smoke/ess/result.rds'),str(ROOT/'reproduction/config/xenium_five_targets.tsv'),str(smoke_out),'1850']
    (smoke_out/'command.json').write_text(json.dumps(smoke_cmd,indent=2)+'\n')
    with (out/'smoke/comparison_retry.log' if a.resume_after_smoke else out/'smoke/comparison.log').open('w') as f:
        smoke_rc=subprocess.run(smoke_cmd,cwd=ROOT,env=env_for(),stdout=f,stderr=subprocess.STDOUT).returncode
    (out/'smoke_status.json').write_text(json.dumps({'status':'PASSED' if smoke_rc==0 else 'FAILED','comparison_exit_code':smoke_rc,
      'interpretation':'Paired small-input smoke only; not Xenium evidence.','comparison_dir':str(smoke_out.relative_to(out))},indent=2)+'\n')
    if smoke_rc!=0:return smoke_rc
    print('Paired smoke passed for both workers; outputs use the common reader/writer and comparator.',flush=True)
    if a.smoke_only:return 0
    # A single active inference process at a time; Original first, then final clean ESS.
    states=[]
    for method in ('original','ess'):
        if available_gib()<MIN_START_GIB:
            state={'method':method,'status':'BLOCKED_RESOURCE_BEFORE_START','MemAvailable_GiB':available_gib()}
            (out/method/'run_status.json').parent.mkdir(parents=True,exist_ok=True)
            (out/method/'run_status.json').write_text(json.dumps(state,indent=2)+'\n')
            states.append(state);break
        state=run_worker(method,INPUT,EXPECTED_INPUT_SHA,out/method,True,a.abort_floor_gib);states.append(state)
        snapshot(out,f'after_{method}')
        if method=='original' and state['status']!='COMPLETED' and available_gib()<MIN_START_GIB:
            blocked={'method':'ess','status':'BLOCKED_RESOURCE_AFTER_original_failure','MemAvailable_GiB':available_gib()}
            (out/'ess').mkdir(exist_ok=True);(out/'ess/run_status.json').write_text(json.dumps(blocked,indent=2)+'\n');states.append(blocked);break
        if state['status'] in ('ABORTED_RESOURCE_GUARD','INTERRUPTED'):break
    (out/'method_run_status.json').write_text(json.dumps(states,indent=2)+'\n')
    if len(states)==2 and all(x.get('status')=='COMPLETED' for x in states):
        if sha(INPUT)!=EXPECTED_INPUT_SHA:raise SystemExit('Input changed after runs; will not compare')
        rc=run_compare(out,627)
        if rc==0:print('Strict full 627-LR comparison completed.',flush=True)
        else:print('Full comparison failed; see comparison.log.',flush=True)
        return rc
    print('One or both full runs did not complete; no paired speed ratio or full agreement claim generated.',flush=True)
    return 2
if __name__=='__main__':raise SystemExit(main())
