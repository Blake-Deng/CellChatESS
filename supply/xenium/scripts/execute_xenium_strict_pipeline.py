#!/usr/bin/env python3
"""Run the strict pair after smoke, then comparison/finalization; retain failures."""
from pathlib import Path
import argparse,csv,hashlib,json,os,signal,subprocess,sys,time,traceback
ROOT=Path(__file__).resolve().parents[1]
def dump(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    out=a.output.resolve();state=out/'pipeline_status.json'
    if state.exists():raise SystemExit('Refusing to restart an existing pipeline')
    dump(out/'pipeline_process.json',{'pid':os.getpid(),'started_epoch':time.time()})
    dump(state,{'status':'RUNNING','stage':'full_sequential_pair','started_epoch':time.time()})
    cmd=[sys.executable,str(ROOT/'reproduction/run_xenium_paired_full627.py'),'--resume-after-smoke','--output',str(out)]
    child=subprocess.Popen(cmd,cwd=ROOT)
    def stop(signum,frame):
        if child.poll() is None:child.send_signal(signum)
    signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
    rc=child.wait()
    dump(out/'pair_driver_exit.json',{'exit_code':rc,'finished_epoch':time.time()})
    if rc:
        dump(state,{'status':'INCOMPLETE','stage':'pair_or_comparison','exit_code':rc,'finished_epoch':time.time(),'full627_pass_claim':False})
        (out/'INCOMPLETE_CN.md').write_text('本次全量配对或比较未完成。请查看各方法 run_status.json、run.log、GNU time 及 comparison.log；不能据此报告完整配对加速比或全量一致。\n')
        return rc
    summary={r['check']:r['value'] for r in csv.DictReader((out/'comparison/full627_comparison_summary.tsv').open(),delimiter='\t')}
    dump(state,{'status':'RUNNING','stage':'finalization','pair_and_chunked_comparison_complete':True})
    cmd=[sys.executable,str(ROOT/'reproduction/finalize_xenium_paired.py'),'--output',str(out)]
    with (out/'finalization.log').open('w') as f:rc=subprocess.run(cmd,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT).returncode
    if rc:
        dump(state,{'status':'INCOMPLETE','stage':'finalization','exit_code':rc,'finished_epoch':time.time()});return rc
    diffs=int(summary['pvalue_different_elements'])
    dump(state,{'status':'RESULTS_READY' if not diffs else 'DIAGNOSIS_REQUIRED','stage':'finished','pvalue_different_elements':diffs,
      'aggregate_counts_identical':summary['aggregate_counts_identical'],'finished_epoch':time.time(),
      'scope':'Only this frozen Xenium input, sources and parameters; no general floating-point identity claim.'})
    return 0
if __name__=='__main__':
    try:sys.exit(main())
    except Exception:
        traceback.print_exc();sys.exit(1)
