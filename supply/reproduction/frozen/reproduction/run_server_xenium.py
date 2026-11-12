#!/usr/bin/env python3
"""Run one full 627-LR PROFILE derivative, tracking the actual R child PID."""
from pathlib import Path
import json,os,subprocess,time,shutil,sys
from server_tools import ROOT,sha,runtime,env1,rows,table,resources
rt=runtime();out=ROOT/'results/new/xenium_profile'
if out.exists():raise SystemExit('Fresh Xenium profile output required')
out.mkdir(parents=True)
mem={l.split(':')[0]:int(l.split()[1])*1024 for l in Path('/proc/meminfo').read_text().splitlines()}
free=shutil.disk_usage(ROOT).free
pre={'available_memory_bytes':mem['MemAvailable'],'disk_free_bytes':free,'required_memory_bytes':180*1024**3,'required_disk_bytes':200*1024**3,'sampling_interval_sec':1.0,'other_campaign_timing_queues':'Must have completed before this stage; see campaign status and process snapshot'}
(out/'resources_before.json').write_text(json.dumps(pre,indent=2)+'\n')
(out/'processes_before.txt').write_text(subprocess.check_output(['ps','-eo','pid,ppid,etime,pcpu,rss,comm'],text=True))
if mem['MemAvailable']<180*1024**3 or free<200*1024**3:
    (out/'completion.json').write_text(json.dumps({'status':'BLOCKED','reason':'Insufficient free memory/storage','resources':pre},indent=2)+'\n');sys.exit(2)
inp=Path('/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/prepared_xenium_ovary_full_grid50.rds');h=sha(inp)
if h!='cc6155a5de127ba34b8a26e7f69e7bffc704448f976107253711623d2a88de37':raise SystemExit('Xenium SHA mismatch '+h)
stage=out/'stage_events.tsv';result=out/'result.rds';lib=ROOT/'server_variants/PROFILE/library'
e=env1(rt);e.update(CELLCHAT_ESS_R_LIB=str(lib),CELLCHAT_COMPAT_R_LIB=rt['cellchat_22_library'],CELLCHATESS_STAGE_FILE=str(stage))
cmd=['/usr/bin/time','-v','-o',str(out/'gnu_time.txt'),rt['rscript'],'--vanilla',str(ROOT/'reproduction/run_xenium_profiled.R'),str(inp),str(result),str(out/'summary.tsv'),str(out/'sessionInfo.txt')]
(out/'command.json').write_text(json.dumps(cmd,indent=2)+'\n');(out/'identity.json').write_text(json.dumps({'input':str(inp),'input_sha256':h,'profile_library':str(lib),'profile_dll_sha256':sha(lib/'CellChatESS/libs/CellChatESS.so'),'source_diff':'server_variants/PROFILE/controlled_diff.patch','candidate_commit':'139a61c33d4c09568f48f4afd0a16dc57bb22f4f'},indent=2)+'\n')
with (out/'run.log').open('w') as f:
    proc=subprocess.Popen(cmd,cwd=ROOT,env=e,stdout=f,stderr=subprocess.STDOUT)
    rpid=None
    while proc.poll() is None and rpid is None:
        pending=[proc.pid];visited=set()
        while pending:
            q=pending.pop()
            if q in visited:continue
            visited.add(q);p=Path('/proc')/str(q)
            try:
                command=p.joinpath('cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace');comm=p.joinpath('comm').read_text().strip()
                if 'run_xenium_profiled.R' in command and comm in ['R','Rscript']:rpid=q;break
                pending.extend(map(int,p.joinpath('task',str(q),'children').read_text().split()))
            except OSError:pass
        if rpid is None:time.sleep(.02)
    sampler=None
    if rpid is not None:
        sampler=subprocess.Popen([sys.executable,str(ROOT/'reproduction/sample_process_memory_server.py'),'--pid',str(rpid),'--output',str(out/'rss_samples.tsv'),'--stages',str(stage),'--interval','1'],cwd=ROOT,env=e,stdout=(out/'sampler.log').open('w'),stderr=subprocess.STDOUT)
        (out/'sampling_pid.json').write_text(json.dumps({'gnu_time_pid':proc.pid,'actual_R_pid':rpid,'start_epoch':time.time()},indent=2)+'\n')
    rc=proc.wait()
    src=sampler.wait(timeout=60) if sampler else -1
complete=rc==0 and src==0
meta={'status':'COMPLETED' if complete else 'FAILED','R_exit_code':rc,'sampler_exit_code':src,'actual_R_pid':rpid,'scope':'ESS full 627-LR absolute resource profile only; no new Original speedup'}
if (out/'gnu_time.txt').exists():meta.update(resources(out/'gnu_time.txt'))
if (out/'summary.tsv').exists():
    meta.update(rows(out/'summary.tsv')[0]);complete=complete and int(meta['n_lr'])==627 and int(meta['cells'])==1156091 and int(meta['groups'])==2357
    meta['status']='COMPLETED' if complete else 'FAILED'
if result.exists():meta['output_path']=str(result);meta['output_file_sha256']=sha(result)
(out/'completion.json').write_text(json.dumps(meta,indent=2)+'\n')
if not complete:sys.exit(1)
print('Xenium full profile completed.',flush=True)
