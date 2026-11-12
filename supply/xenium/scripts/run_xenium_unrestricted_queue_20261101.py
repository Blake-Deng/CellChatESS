#!/usr/bin/env python3
"""Run a fresh strict Xenium Original/ESS pair without application-side memory thresholds."""
from pathlib import Path
import argparse, csv, hashlib, json, os, shutil, subprocess, sys, time
ROOT=Path(__file__).resolve().parents[1]
SOURCE_RUN=ROOT/'results/new/xenium_paired_full627_20261101_strict03'
RUNNER=ROOT/'reproduction/run_xenium_paired_full627_unrestricted_20261101.py'
def available_gib():
    line=next(x for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:'))
    return int(line.split()[1])/1024**2
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()
def write_json(p,x):
    t=p.with_suffix('.tmp');t.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n');t.replace(p)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--control',type=Path,required=True);ap.add_argument('--attempt',type=Path,required=True);a=ap.parse_args()
    control=a.control.resolve();out=a.attempt.resolve();control.mkdir(parents=True,exist_ok=True)
    if out.exists(): raise SystemExit('Refusing to overwrite '+str(out))
    lock=control/'.queue.lock'
    fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600);os.write(fd,str(os.getpid()).encode());os.close(fd)
    start=time.time();mem0=available_gib()
    with (control/'queue.log').open('a',buffering=1) as log:
        log.write(f'Queue started {time.strftime("%Y-%m-%d %H:%M:%S %Z")} with MemAvailable={mem0:.1f} GiB; no application memory start/abort thresholds.\n')
        out.mkdir(parents=True)
        for d in ('sources','environment','audit','auxiliary','scripts'):
            shutil.copytree(SOURCE_RUN/d,out/d)
        (out/'smoke').mkdir()
        for method in ('original','ess'): shutil.copytree(SOURCE_RUN/'smoke'/method,out/'smoke'/method)
        shutil.copy2(SOURCE_RUN/'smoke_status.json',out/'smoke_status.json')
        shutil.copy2(ROOT/'reproduction/audit_xenium_paired_delivery.py',out/'scripts/audit_xenium_paired_delivery.py')
        shutil.copy2(RUNNER,out/'scripts/run_xenium_paired_full627_unrestricted_20261101.py')
        shutil.copy2(Path(__file__),out/'scripts/run_xenium_unrestricted_queue_20261101.py')
        (out/'RETRY_CONTEXT_CN.md').write_text(
            '# 无应用层内存阈值的严格重跑\n\n'
            '此前 strict03 因自定义 MemAvailable 180 GiB 中止阈值停止：Original 进度条到 100%，函数未返回且未保存完整结果；ESS 未运行。该旧尝试保持原样。\n\n'
            '本目录由作者明确要求启动的新配对重跑产生。入口副本仅将应用层内存启动阈值和中止阈值设为 0；未修改 CellChat/ESS 源码、输入、参数、线程数、计时或统计规则。系统/操作系统仍可能因物理内存不足终止进程。\n\n'
            '复用的只有同一冻结 worker 已通过的 PBMC smoke 输出；Xenium Original 和 ESS 均从头、独立进程、顺序运行。\n')
        write_json(control/'queue_status.json',{'status':'RUNNING_STRICT_PAIR','queue_pid':os.getpid(),'driver_pid':None,
            'attempt_output':str(out),'started_epoch':start,'MemAvailable_GiB_at_launch':mem0,
            'application_memory_start_threshold_GiB':0,'application_memory_abort_threshold_GiB':0,
            'input_sha256':'cc6155a5de127ba34b8a26e7f69e7bffc704448f976107253711623d2a88de37',
            'source_script_sha256':sha(RUNNER),'source_script_diff':str(control/'runner_memory_thresholds.patch')})
        cmd=[sys.executable,str(RUNNER),'--resume-after-smoke','--output',str(out)]
        (control/'run_command.json').write_text(json.dumps(cmd,indent=2)+'\n')
        with (control/'driver.log').open('w') as driver_log:
            child=subprocess.Popen(cmd,cwd=ROOT,stdout=driver_log,stderr=subprocess.STDOUT,start_new_session=True)
            write_json(control/'driver_process.json',{'pid':child.pid,'command':cmd,'started_epoch':time.time()})
            last=0
            while child.poll() is None:
                now=time.time()
                if now-last>=60:
                    state={'status':'RUNNING_STRICT_PAIR','queue_pid':os.getpid(),'driver_pid':child.pid,
                           'attempt_output':str(out),'updated_epoch':now,'elapsed_sec':now-start,
                           'MemAvailable_GiB':available_gib(),'no_application_memory_thresholds':True}
                    write_json(control/'queue_status.json',state)
                    log.write(f'{time.strftime("%Y-%m-%d %H:%M:%S %Z")} elapsed={now-start:.0f}s MemAvailable={state["MemAvailable_GiB"]:.1f} GiB driver_pid={child.pid}\n');last=now
                time.sleep(15)
            rc=child.wait()
        status_file=out/'method_run_status.json'
        methods=json.loads(status_file.read_text()) if status_file.is_file() else []
        completed=(rc==0 and len(methods)==2 and all(x.get('status')=='COMPLETED' for x in methods)
                   and (out/'comparison_status.json').is_file()
                   and json.loads((out/'comparison_status.json').read_text()).get('status')=='COMPLETED')
        if completed:
            write_json(out/'pipeline_status.json',{'status':'RUNNING','stage':'finalization','pair_comparison_complete':True,'updated_epoch':time.time()})
            with (out/'finalization.log').open('w') as f:
                fin=subprocess.run([sys.executable,str(ROOT/'reproduction/finalize_xenium_paired.py'),'--output',str(out)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
            completed=fin.returncode==0
            summary={r['check']:r['value'] for r in csv.DictReader((out/'comparison/full627_comparison_summary.tsv').open(),delimiter='\t')}
            state='RESULTS_READY' if completed and int(summary['pvalue_different_elements'])==0 else ('DIAGNOSIS_REQUIRED' if completed else 'INCOMPLETE_FINALIZATION')
            write_json(out/'pipeline_status.json',{'status':state,'stage':'finished','pvalue_different_elements':summary.get('pvalue_different_elements'),
                'aggregate_counts_identical':summary.get('aggregate_counts_identical'),'finished_epoch':time.time(),
                'scope':'Only this frozen Xenium input, sources, parameters and actual conditions.'})
        else:
            write_json(out/'pipeline_status.json',{'status':'INCOMPLETE','stage':'pair_or_comparison','driver_exit_code':rc,'finished_epoch':time.time()})
            (out/'INCOMPLETE_CN.md').write_text('本次全量配对或比较未完成。原始日志、GNU time、RSS 和阶段记录均保留；不得据此声称全量验证通过或报告完整配对加速比。\n')
        audit=subprocess.run([sys.executable,str(ROOT/'reproduction/audit_xenium_paired_delivery.py'),'--output',str(out)],cwd=ROOT,
                             stdout=(control/'terminal_audit.log').open('w'),stderr=subprocess.STDOUT)
        terminal='COMPLETED_EVIDENCE_INDEXED' if completed and audit.returncode==0 else 'INCOMPLETE_EVIDENCE_INDEXED'
        write_json(control/'queue_status.json',{'status':terminal,'queue_pid':os.getpid(),'driver_pid':child.pid,'attempt_output':str(out),
            'driver_exit_code':rc,'pair_and_comparison_completed':completed,'terminal_audit_exit_code':audit.returncode,
            'finished_epoch':time.time(),'elapsed_sec':time.time()-start,'MemAvailable_GiB_at_end':available_gib(),
            'no_application_memory_thresholds':True})
        log.write(f'Queue ended {terminal}; pair_rc={rc}, audit_rc={audit.returncode}.\n')
        return 0 if completed and audit.returncode==0 else 2
if __name__=='__main__': raise SystemExit(main())
