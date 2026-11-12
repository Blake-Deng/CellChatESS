#!/usr/bin/env python3
"""Correct the input selection explicitly; preserve earlier high-retention records."""
from server_tools import *
rt=runtime();out=ROOT/'results/new/highest_retention_CC22'
if out.exists():raise SystemExit('Choose fresh addendum output')
out.mkdir(parents=True)
f=rows(ROOT/'results/new/core22_statistics/dataset_characteristics_216.tsv')
x=max(f,key=lambda r:float(r['f_LR']))
assert x['experiment_id']=='3CA_data__SmartSeq2__cells-1000__rep-2__engine-both__ablation-full'
assert sha(x['input_file'])==x['input_sha256']
table(out/'selected_input.tsv',[x])
libs={'A0':Path(rt['ess_library']),'A4':ROOT/'server_variants/A4_cache_no_pruning/library','A5':Path(rt['ess_library'])}
results=[];comparisons=[]
def run(stage,rep,phase,pos):
 d=out/phase/f'rep{rep}'/stage
 rc=record_run([rt['rscript'],'--vanilla',ROOT/'reproduction/run_ablation_worker.R',stage,x['input_file'],d/'output',libs[stage],rt['cellchat_22_library']],d,env1(rt))
 if rc:raise RuntimeError(str(d))
 r=rows(d/'output/measurement.tsv')[0];r.update(experiment_id=x['experiment_id'],dataset_id=x['dataset_id'],scale=x['scale'],replicate=rep,stage=stage,phase=phase,order_position=pos,input_sha256=x['input_sha256'],f_LR=x['f_LR'],exit_code=rc)
 dll=Path(rt['cellchat_22_library'])/'CellChat/libs/CellChat.so' if stage=='A0' else libs[stage]/'CellChatESS/libs/CellChatESS.so'
 r.update(resources(d/'gnu_time.txt'));r.update(dll_sha256=sha(dll),result_path=str(d/'output/result.rds'));results.append(r);table(out/'all_measurements.tsv',results);print('DONE',stage,rep,phase,flush=True)
 return d/'output/result.rds'
def check(ref,y,stage,rep,phase):
 d=y.parent.parent/'numerical_comparison';ok=compare(ref,y,d,rt)
 for r in rows(d/'comparison.tsv'):
  r.update(experiment_id=x['experiment_id'],stage=stage,replicate=rep,phase=phase,status='PASS' if ok else 'FAIL');comparisons.append(r)
 table(out/'numerical_comparisons.tsv',comparisons)
 if not ok:raise RuntimeError('Numerical comparison failed')
ref=run('A0',0,'validation',1)
for pos,stage in enumerate(['A4','A5'],2):check(ref,run(stage,0,'validation',pos),stage,0,'validation')
for rep,order in enumerate([['A4','A5'],['A5','A4'],['A4','A5']],1):
 for pos,stage in enumerate(order,1):check(ref,run(stage,rep,'timing',pos),stage,rep,'timing')
(out/'completion.json').write_text(json.dumps({'status':'COMPLETED','purpose':'Highest retained-fraction input addendum; earlier rep1 retained under ablation_CC22','f_LR':x['f_LR'],'experiment_id':x['experiment_id'],'validation_processes':3,'timing_processes':6,'CellChat':'2.2.0.9001','ESS':'0.2.2','within_experiment_concurrency':1},indent=2)+'\n')
