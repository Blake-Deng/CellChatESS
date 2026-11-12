#!/usr/bin/env python3
"""Run frozen-source server verification and then isolated experiment phases.

All output paths are new. No historical record or public repository is modified.
Use --prepare-only for build/validation; use --experiments-only after that passes.
"""
from pathlib import Path
import argparse, csv, datetime, hashlib, json, os, re, shutil, subprocess, sys, time

ROOT=Path(__file__).resolve().parents[1]
NEW=ROOT/'results/new'
P=argparse.ArgumentParser()
P.add_argument('--prepare-only',action='store_true')
P.add_argument('--experiments-only',action='store_true')
P.add_argument('--primary-jobs',type=int,default=8)
P.add_argument('--autozyme-jobs',type=int,default=4)
A=P.parse_args()
rt=json.loads((ROOT/'server_install/runtime.json').read_text())
R=Path(rt['rscript']).with_name('R');RS=rt['rscript']
env=os.environ.copy()
env['PATH']=str(R.parent)+os.pathsep+env.get('PATH','')
env['R_LIBS']=rt['dependency_library']
env['AUTOZYME_DISABLED']='true'
env['MAKEFLAGS']='-j2'
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS',
          'VECLIB_MAXIMUM_THREADS','RCPP_PARALLEL_NUM_THREADS','AUTOZYME_THREADS','ZYME_THREADS']:
    env[k]='1'
env['OMP_DYNAMIC']='FALSE'
env['OMP_THREAD_LIMIT']='4'
state_path=ROOT/'server_campaign_status.json'
state=json.loads(state_path.read_text()) if state_path.exists() else {'tasks':{},'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
state.update(pid=os.getpid(),baseline='CellChat 2.2.0.9001',ESS_commit='139a61c33d4c09568f48f4afd0a16dc57bb22f4f',
             primary_jobs=A.primary_jobs,autozyme_jobs=A.autozyme_jobs,publications='NOT_AUTHORIZED')

def save():
    state['updated_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    tmp=state_path.with_suffix('.tmp');tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(state_path)

def run(name,command,*,extra_env=None,cwd=ROOT,required=False):
    previous=state['tasks'].get(name,{})
    if previous.get('status')=='PASS':
        print('ALREADY PASS',name,flush=True);return True
    log=ROOT/'logs'/(name+'.log')
    if log.exists():
        log=ROOT/'logs'/(name+f'_attempt_{int(time.time())}.log')
    ee=env.copy();ee.update(extra_env or {})
    entry={'status':'RUNNING','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
           'command':[str(z) for z in command],'cwd':str(cwd),'log':str(log.relative_to(ROOT))}
    state['tasks'][name]=entry;state['current_task']=name;save()
    print('START',name,flush=True)
    with log.open('w') as f:
        proc=subprocess.Popen([str(z) for z in command],cwd=cwd,env=ee,stdout=f,stderr=subprocess.STDOUT)
        entry['child_pid']=proc.pid;save();rc=proc.wait()
    entry.update(exit_code=rc,status='PASS' if rc==0 else 'FAIL',finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());save()
    print(entry['status'],name,flush=True)
    if rc and required:raise RuntimeError(name+' failed: '+str(log))
    return rc==0

def fresh(p):
    if p.exists():raise RuntimeError('Fresh path required: '+str(p))
    p.mkdir(parents=True)

def preparation():
    # These were executed before this controller was installed. Verify their outputs.
    pre=NEW/'version_audit/PASS.txt'
    smoke=NEW/'server_pbmc_smoke/checks.tsv'
    if not pre.exists() or not smoke.exists():raise RuntimeError('Complete preflight and PBMC smoke first')
    checks=list(csv.DictReader(smoke.open(),delimiter='\t'))
    if not checks or any(r['status']!='PASS' for r in checks):raise RuntimeError('PBMC smoke did not pass')
    state['tasks']['source_and_version_preflight']={'status':'PASS','evidence':str(pre.relative_to(ROOT))}
    state['tasks']['three_method_real_smoke']={'status':'PASS','evidence':str(smoke.relative_to(ROOT)),'n_checks':len(checks)};save()
    run('capture_server_environment',[RS,'--vanilla','reproduction/capture_server_environment.R',
        'server_install/runtime.json','environment/server'],required=True)
    checks_dir=ROOT/'server_build_checks';checks_dir.mkdir(exist_ok=True)
    run('build_final_source',[R,'CMD','build','--no-manual','--no-build-vignettes',ROOT/'software/CellChatESS_candidate_git'],cwd=checks_dir,required=True)
    buildenv={'R_LIBS':os.pathsep.join([rt['cellchat_22_library'],rt['ess_library'],rt['dependency_library']]),
              '_R_CHECK_FORCE_SUGGESTS_':'false'}
    run('R_CMD_check',[R,'CMD','check','--no-manual','--no-build-vignettes','CellChatESS_0.2.2.tar.gz'],
        cwd=checks_dir,extra_env=buildenv,required=True)
    run('synthetic_three_method_smoke',[sys.executable,'validation/run_three_method_smoke.py','--runtime',
        'server_install/runtime.json','--input','validation/local_smoke/data/synthetic_structures_CC22.rds',
        '--output',NEW/'synthetic_smoke'],required=True)
    run('CC22_alignment',[RS,'--vanilla','validation/check_CC22_alignment.R',NEW/'CC22_alignment',
        rt['ess_library'],rt['cellchat_22_library'],rt['dependency_library']],required=True)
    real_input=list(csv.DictReader((ROOT/'reproduction/config/primary_216_CC22.tsv').open(),delimiter='\t'))[0]['input_file']
    for version,key in [('2.2.0.9001','cellchat_22_library'),('1.6.1','cellchat_16_library')]:
        tag='cc22' if version=='2.2.0.9001' else 'cc16'
        e={'CELLCHAT_ESS_R_LIB':rt['ess_library'],'CELLCHAT_COMPAT_R_LIB':rt[key],
           'CELLCHAT_TEST_INPUT':real_input,'CELLCHAT_REFERENCE_VERSION':version}
        run('structure_26_'+tag,[RS,'--vanilla','validation/run_structure_tests_server.R',NEW/('structure_'+tag)],extra_env=e,required=True)
        run('input_guards_5_'+tag,[RS,'--vanilla','validation/run_input_guards_server.R',real_input,
            NEW/('input_guards_'+tag+'.tsv')],extra_env=e,required=True)
    # Controlled builds are installed before any full timing queue begins.
    libs={}
    for name in ['A1_dense_full_gene','A3_no_cache_no_pruning','A4_cache_no_pruning','PROFILE','TRACE']:
        lib=ROOT/'server_variants'/name/'library';lib.mkdir(exist_ok=True)
        run('install_variant_'+name,[R,'CMD','INSTALL','--preclean','--clean','--library='+str(lib),
            ROOT/'server_variants'/name/'CellChatESS'],required=True)
        libs[name]=str(lib)
        out=NEW/('variant_smoke_'+name);fresh(out)
        method='dense' if name.startswith('A1_') else 'ess'
        run('smoke_variant_'+name,[RS,'--vanilla','validation/smoke_worker.R',method,
            ROOT/'validation/local_smoke/data/synthetic_structures_CC22.rds',out,lib,
            rt['cellchat_22_library'],rt['autozyme_library'],rt['dependency_library'],'2.2.0.9001','TRUE'],required=True)
        run('compare_variant_'+name,[RS,'--vanilla','reproduction/compare_results.R',
            NEW/'synthetic_smoke/original/result.rds',out/'result.rds',out/'comparison.tsv'],required=True)
    (ROOT/'server_variants/runtime.json').write_text(json.dumps(libs,indent=2)+'\n')
    probe=ROOT/'environment/server/openmp_team_probe.so'
    run('compile_openmp_team_probe',[str(R.parent/'x86_64-conda-linux-gnu-cc'),'-shared','-fPIC','-O2','-fopenmp',
        ROOT/'reproduction/openmp_team_probe.c','-ldl','-o',probe],required=True)
    for version,key in [('2.2.0.9001','cellchat_22_library'),('1.6.1','cellchat_16_library')]:
        tag='cc22' if version=='2.2.0.9001' else 'cc16'
        root=NEW/('openmp_'+tag);fresh(root)
        team=root/'actual_openmp_teams.tsv'
        team.write_text('pid\tthreads_requested\tactual_team_size\tcaller\n')
        for n in [1,4]:
            run('openmp_'+tag+'_'+str(n),[RS,'--vanilla','validation/check_openmp_server.R',real_input,
                root/('threads'+str(n)),rt['ess_library'],rt[key],version,str(n)],
                extra_env={'LD_PRELOAD':str(probe),'CELLCHATESS_OPENMP_TEAM_LOG':str(team)},required=True)
        run('openmp_compare_'+tag,[RS,'--vanilla','reproduction/compare_results.R',
            root/'threads1/result.rds',root/'threads4/result.rds',root/'comparison.tsv'],required=True)
        teams=list(csv.DictReader(team.open(),delimiter='\t'))
        own=[r for r in teams if 'CellChatESS.so' in r['caller']]
        if not own or any(int(r['actual_team_size'])!=4 for r in own):
            state['tasks']['openmp_actual_team_'+tag]={'status':'FAIL','evidence':str(team.relative_to(ROOT))};save()
            raise RuntimeError('Failed to directly verify real OpenMP team size 4 in ESS')
        state['tasks']['openmp_actual_team_'+tag]={'status':'PASS','teams_observed':len(own),
            'actual_team_sizes':sorted({int(r['actual_team_size']) for r in own}),'evidence':str(team.relative_to(ROOT))};save()
    run('primary216_worker_smoke',[sys.executable,'reproduction/run_unified.py','primary216',
        '--runtime','server_install/runtime.json','--output',NEW/'primary_smoke','--limit','1','--jobs','1','--execute'],required=True)
    state['preparation']='PASS';state['current_task']='PREPARATION_COMPLETE';save()

def experiments():
    if state.get('preparation')!='PASS':raise RuntimeError('Preparation gates have not all passed')
    # No new builds or other in-house CPU-heavy analyses run concurrently with these queues.
    run('primary216_full',[sys.executable,'-u','reproduction/run_unified.py','primary216',
        '--runtime','server_install/runtime.json','--output',NEW/'primary216','--jobs',str(A.primary_jobs),'--execute'])
    run('autozyme192_worker_smoke',[sys.executable,'-u','reproduction/run_unified.py','autozyme192',
        '--runtime','server_install/runtime.json','--output',NEW/'autozyme_smoke','--limit','1','--jobs','1','--execute'],required=True)
    run('autozyme192_full',[sys.executable,'-u','reproduction/run_unified.py','autozyme192',
        '--runtime','server_install/runtime.json','--output',NEW/'autozyme192','--jobs',str(A.autozyme_jobs),'--execute'])
    # Independent follow-ups continue even if an individual queue records numerical differences.
    for name,args in [
        ('ablation','run_server_ablation.py'),('xenium_target_traces','run_server_traces.py'),
        ('xenium_full_profile','run_server_xenium.py'),('final_statistics_figures','finalize_server_results.py')]:
        script=ROOT/'reproduction'/args
        if script.exists():run(name,[sys.executable,'-u',script])
        else:state['tasks'][name]={'status':'NOT_RUN','reason':'Server follow-up entrypoint not yet prepared','expected_script':str(script)};save()
    state['current_task']='FINISHED';state['status']='COMPLETE' if all(r['status']=='PASS' for r in state['tasks'].values()) else 'COMPLETE_WITH_UNRESOLVED_ITEMS';save()

try:
    if not A.experiments_only:preparation()
    if not A.prepare_only:experiments()
except BaseException as e:
    state['status']='STOPPED_AT_VALIDATION_GATE';state['error']=repr(e);save();raise
