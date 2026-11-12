#!/usr/bin/env python3
"""Diagnostic traces for actual full-run p-value differences, after timing ends."""
from pathlib import Path
import argparse,csv,hashlib,json,os,subprocess,tempfile,shutil,time
ROOT=Path(__file__).resolve().parents[1]
TRACE_SHA='713a1a47eab864ae8e01f9c3cb4331003ee5fad5c04c3bcee54d5ae8a2465750'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for x in iter(lambda:f.read(8*1024*1024),b''):h.update(x)
    return h.hexdigest()
def table(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();out=a.output.resolve()
    differences=list(csv.DictReader((out/'comparison/pvalue_differences.tsv').open(),delimiter='\t'))
    dest=out/'difference_diagnostics'
    if dest.exists():raise SystemExit('Refusing to overwrite diagnosis')
    dest.mkdir()
    if not differences:
        (dest/'status.json').write_text(json.dumps({'status':'NOT_REQUIRED','reason':'No full-run p-value differences; auxiliary historical five-target traces are separate evidence.'},indent=2)+'\n');return
    rt=json.loads((ROOT/'server_install/runtime.json').read_text())
    trace=ROOT/'server_variants/TRACE';dll=trace/'library/CellChatESS/libs/CellChatESS.so'
    if sha(dll)!=TRACE_SHA:raise SystemExit('Trace DLL identity mismatch')
    patch=trace/'controlled_diff.patch'
    with tempfile.TemporaryDirectory(prefix='xenium_trace_identity_') as t:
        check=Path(t)/'source';shutil.copytree(ROOT/'software/CellChatESS',check)
        run=subprocess.run(['patch','--batch','-p1','-i',str(patch)],cwd=check,capture_output=True,text=True)
        (dest/'patch_reapplication.log').write_text(run.stdout+run.stderr)
        if run.returncode:raise SystemExit('Trace patch cannot be reapplied to final source')
        for rel in ['DESCRIPTION','R/accel_cellchat.R','R/utils.R','src/cellchat_prob_cpp.cpp','src/RcppExports.cpp']:
            if sha(check/rel)!=sha(trace/'CellChatESS'/rel):raise SystemExit('Trace source differs from declared pure patch: '+rel)
    (dest/'source_identity.json').write_text(json.dumps({'candidate_commit':'139a61c33d4c09568f48f4afd0a16dc57bb22f4f','trace_DLL_sha256':sha(dll),'pure_trace_patch_sha256':sha(patch),'patch_reapplied_and_R_cpp_verified':True,'timing_use':'diagnostic only; excluded from paired performance measurements'},indent=2)+'\n')
    env=os.environ.copy();env.update({'R_LIBS':rt['dependency_library'],'R_LIBS_USER':'','R_LIBS_SITE':'','R_PROFILE_USER':'/dev/null','R_ENVIRON_USER':'/dev/null','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','BLIS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','AUTOZYME_DISABLED':'true','CELLCHAT_REFERENCE_VERSION':'2.2.0.9001','LC_ALL':'C'})
    original_meta=json.loads((out/'original/run_metadata.json').read_text());input_path=Path(original_meta['input']['input_path'])
    if sha(input_path)!=original_meta['input']['input_sha256']:raise SystemExit('Input hash changed before diagnosis')
    results=[]
    for index,lr in enumerate(dict.fromkeys(r['lr_name'] for r in differences),1):
        selected=[r for r in differences if r['lr_name']==lr]
        batch=dest/f'LR_{index:04d}';batch.mkdir()
        targets=[{k:r[k] for k in ('lr_name','source_group','target_group')} for r in selected]
        table(batch/'targets.tsv',targets)
        cmd=[rt['rscript'],'--vanilla',str(ROOT/'reproduction/run_target_trace_CC22.R'),str(input_path),str(batch/'targets.tsv'),str(batch/'trace'),str(trace/'library'),rt['cellchat_22_library']]
        (batch/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
        with (batch/'run.log').open('w') as f:rc=subprocess.run(cmd,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
        result={'lr_name':lr,'targets':len(targets),'exit_code':rc,'status':'COMPLETED' if rc==0 else 'FAILED','evidence':str(batch.relative_to(out)),'diagnostic_scope':'one LR, all cells/groups, global normalization anchor, 100 bootstraps per target; not a repeat of the full performance run'}
        if rc==0:
            observed=list(csv.DictReader((batch/'trace/target_summary.tsv').open(),delimiter='\t'))
            lookup={(r['lr_name'],r['source_group'],r['target_group']):r for r in observed}
            matches=[]
            for r in selected:
                t=lookup[(r['lr_name'],r['source_group'],r['target_group'])]
                matches.append({'lr_name':lr,'source_group':r['source_group'],'target_group':r['target_group'],
                  'original_p_matches_full_run':float(t['original_pvalue'])==float(r['original_value']),
                  'ESS_p_matches_full_run':float(t['ESS_pvalue'])==float(r['ESS_value']),
                  'original_observed_probability_hex_matches_full_run':t['original_Pobs_hex']==r['original_observed_probability_hex'],
                  'ESS_observed_probability_hex_matches_full_run':t['ESS_Pobs_hex']==r['ESS_observed_probability_hex']})
            table(batch/'trace_vs_full_run.tsv',matches)
            result['all_target_trace_values_match_full_run']=all(all(v for k,v in z.items() if k.endswith('matches_full_run')) for z in matches)
        results.append(result)
        (dest/'status.json').write_text(json.dumps({'status':'RUNNING','results':results,'updated_epoch':time.time()},indent=2)+'\n')
    (dest/'status.json').write_text(json.dumps({'status':'COMPLETED' if all(r['status']=='COMPLETED' and r.get('all_target_trace_values_match_full_run') for r in results) else 'INCOMPLETE','results':results,'note':'Strict comparison and probability formulas were not altered; causes require review of the saved direct traces.'},indent=2)+'\n')
if __name__=='__main__':main()
