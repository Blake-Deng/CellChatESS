#!/usr/bin/env python3
"""Five real Xenium targets, first CC16 then CC22, final TRACE derivative."""
import json,sys
from server_tools import ROOT,sha,rows,table,runtime,env1,record_run
rt=runtime();out=ROOT/'results/new/xenium_target_trace'
if out.exists():raise SystemExit('Fresh trace output required')
out.mkdir(parents=True)
inp='/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/prepared_xenium_ovary_full_grid50.rds'
expected='cc6155a5de127ba34b8a26e7f69e7bffc704448f976107253711623d2a88de37'
actual=sha(inp)
if actual!=expected:raise SystemExit('Xenium SHA256 mismatch '+actual)
lib=ROOT/'server_variants/TRACE/library';records=[]
for version,key,tag in [('1.6.1','cellchat_16_library','cc16'),('2.2.0.9001','cellchat_22_library','cc22')]:
    d=out/tag;e=env1(rt);e['CELLCHAT_REFERENCE_VERSION']=version
    cmd=[rt['rscript'],'--vanilla',ROOT/'reproduction/run_target_trace_CC22.R',inp,ROOT/'reproduction/config/xenium_five_targets.tsv',d/'output',lib,rt[key]]
    rc=record_run(cmd,d,e)
    row={'CellChat':version,'ESS':'0.2.2','build':'TRACE','input_sha256':actual,'dll_sha256':sha(lib/'CellChatESS/libs/CellChatESS.so'),'exit_code':rc,'status':'FAILED','evidence':str(d.relative_to(ROOT))}
    status=d/'output/trace_status.tsv'
    if status.exists():
        row.update(rows(status)[0])
        row['status']=('COMPLETED' if rc==0 and int(row['comparisons'])==500 and int(row['different_decisions'])==0 and int(row['pvalue_different_targets'])==0 and row['permutations_equal']=='TRUE'
                       else ('COMPLETED_WITH_DIFFERENCES' if rc==0 and int(row['comparisons'])==500 else 'FAILED'))
    records.append(row);table(out/'validation_summary.tsv',records)
    print(tag,row,flush=True)
(out/'completion.json').write_text(json.dumps({'status':('COMPLETED' if all(r['status']=='COMPLETED' for r in records) else ('COMPLETED_WITH_DIFFERENCES' if all(r['status'] in ('COMPLETED','COMPLETED_WITH_DIFFERENCES') for r in records) else 'FAILED')),'records':records,'scope':'Five actual Xenium entries, four LR rows, 100 bootstrap decisions each per CellChat version. Not a 627-LR numerical comparison.'},indent=2)+'\n')
if any(r['status']=='FAILED' for r in records):sys.exit(1)
