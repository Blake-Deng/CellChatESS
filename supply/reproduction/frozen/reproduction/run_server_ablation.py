#!/usr/bin/env python3
"""Controlled CC22 ablation: full-input validation, then three timing repeats."""
from pathlib import Path
import json, sys, time
from server_tools import ROOT,sha,rows,table,runtime,env1,resources,record_run,compare
rt=runtime();out=ROOT/'results/new/ablation_CC22'
if out.exists():raise SystemExit('Fresh ablation output required')
out.mkdir(parents=True)
manifest=rows(ROOT/'reproduction/config/ablation_inputs_CC22.tsv')
libs={s:Path(rt['ess_library']) for s in ['A0','A1','A2','A3','A4','A5']}
libs.update(A1=ROOT/'server_variants/A1_dense_full_gene/library',A3=ROOT/'server_variants/A3_no_cache_no_pruning/library',A4=ROOT/'server_variants/A4_cache_no_pruning/library')
results=[];comparisons=[];failures=[]
for x in manifest:
    if sha(x['input_file'])!=x['input_sha256']:raise SystemExit('Input hash mismatch '+x['experiment_id'])
table(out/'input_hash_verification.tsv',[dict(experiment_id=x['experiment_id'],input_file=x['input_file'],sha256=x['input_sha256'],status='PASS') for x in manifest])

def run_stage(x,stage,rep,phase,position):
    d=out/x['experiment_id']/phase/f'rep{rep}'/stage
    cmd=[rt['rscript'],'--vanilla',ROOT/'reproduction/run_ablation_worker.R',stage,x['input_file'],d/'output',libs[stage],rt['cellchat_22_library']]
    rc=record_run(cmd,d,env1(rt));meta=dict(experiment_id=x['experiment_id'],dataset_id=x['dataset_id'],scale=x['scale'],replicate=rep,stage=stage,phase=phase,order_position=position,input_sha256=x['input_sha256'],exit_code=rc)
    if rc:
        failures.append({**meta,'evidence':str(d.relative_to(ROOT))});return None
    rec=rows(d/'output/measurement.tsv')[0]
    dll=(Path(rt['cellchat_22_library'])/'CellChat/libs/CellChat.so') if stage=='A0' else libs[stage]/'CellChatESS/libs/CellChatESS.so'
    rec.update(meta);rec.update(resources(d/'gnu_time.txt'));rec.update(dll_sha256=sha(dll),result_path=str(d/'output/result.rds'))
    results.append(rec);table(out/'all_measurements.tsv',results)
    print('DONE',phase,x['experiment_id'],rep,stage,flush=True)
    return d/'output/result.rds'

def check(reference,result,x,stage,rep,phase):
    if result is None:return False
    d=result.parent.parent/'numerical_comparison';ok=compare(reference,result,d,rt)
    if (d/'comparison.tsv').exists():
        for r in rows(d/'comparison.tsv'):
            r.update(experiment_id=x['experiment_id'],phase=phase,replicate=rep,stage=stage,status='PASS' if ok else 'FAIL');comparisons.append(r)
        table(out/'numerical_comparisons.tsv',comparisons)
    if not ok:failures.append(dict(experiment_id=x['experiment_id'],stage=stage,phase=phase,replicate=rep,reason='numerical comparison failed',evidence=str(d.relative_to(ROOT))))
    return ok

stages=['A0','A1','A2','A3','A4','A5']
orders=[stages,list(reversed(stages)),['A3','A4','A5','A0','A1','A2']]
for x in [r for r in manifest if r['role']=='A0_A5_representative']:
    validation={s:run_stage(x,s,0,'validation',i+1) for i,s in enumerate(stages)}
    ref=validation['A0']
    ok=ref is not None
    for s in stages[1:]:ok=check(ref,validation[s],x,s,0,'validation') and ok if ref else False
    if not ok:continue
    for rep,order in enumerate(orders,1):
        for pos,s in enumerate(order,1):check(ref,run_stage(x,s,rep,'timing',pos),x,s,rep,'timing')

for x in [r for r in manifest if r['role']=='high_retention_A4_A5']:
    ref=run_stage(x,'A0',0,'high_retention_validation',1)
    valid=True
    for pos,s in enumerate(['A4','A5'],2):
        result=run_stage(x,s,0,'high_retention_validation',pos)
        valid=(check(ref,result,x,s,0,'high_retention_validation') if ref else False) and valid
    if not valid:continue
    for rep,order in enumerate([['A4','A5'],['A5','A4'],['A4','A5']],1):
        for pos,s in enumerate(order,1):check(ref,run_stage(x,s,rep,'high_retention_timing',pos),x,s,rep,'high_retention_timing')
expected={'validation':18,'timing':54,'high_retention_validation':3,'high_retention_timing':6}
actual={p:sum(r['phase']==p for r in results) for p in expected}
complete=not failures and actual==expected
(out/'completion.json').write_text(json.dumps({'status':'COMPLETED' if complete else 'FAILED','expected_method_processes':expected,'actual_method_processes':actual,'failures':failures,'baseline':'2.2.0.9001','candidate_commit':'139a61c33d4c09568f48f4afd0a16dc57bb22f4f','within_experiment_concurrency':1,'validation_excluded_from_performance_statistics':True,'stage_orders':orders},indent=2)+'\n')
if not complete:sys.exit(1)
