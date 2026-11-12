"""Cross-phase report kept with the fresh run, reporting all independent anomaly categories."""
import argparse,json,csv,hashlib
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root
r=root/'results';cross=r/'crosscheck'
def read(path):return json.loads(path.read_text())
def tsv(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
stagefile=cross/'phase_log.json'; phase=read(stagefile)
report={'status':'PASS' if all(e['exit_code']==0 for e in phase) else 'REQUIRES_REVIEW','stages':phase}
for key,path in {'core':cross/'core/audit.json','legacy':r/'legacy_cc16_216/completion.json','preparation':cross/'preparation129/completion.json','xenium_resource_ratios':r/'xenium/fresh_ratios.json','structure_and_threads':r/'auxiliary/validation_status.json','trace':r/'auxiliary/xenium_target_trace/archive_comparison_summary.json','s7_input_features':r/'auxiliary/core22_statistics/comparison_with_archive.json','ablation':r/'auxiliary/ablation_CC22/completion.json','high_retention':r/'auxiliary/highest_retention_CC22/completion.json'}.items():
 if path.is_file():report[key]=read(path)
 else:report[key]={'status':'MISSING','expected':str(path)};report['status']='REQUIRES_REVIEW'
if (r/'figures/statistics.json').is_file():
 st=read(r/'figures/statistics.json');report['new_performance_statistics']={k:st[k] for k in ['primary','autozyme','deduplicated_primary','feature_associations','peak_memory_blocks_ESS_lower','median_peak_RSS_GiB','panel_F_max_absolute_difference','panel_J_points_per_method']}
if (cross/'preparation129/results.tsv').is_file():
 s=tsv(cross/'preparation129/results.tsv');report['preparation_serialization']={'sha_identical':sum(v['expected_sha256']==v['fresh_serialized_sha256'] for v in s),'compared_inputs':len(s),'semantic_agreement':sum(v['group_and_array_slot_agreement']=='True' for v in s)}
if (r/'legacy_cc16_216/results.tsv').is_file():
 z=tsv(r/'legacy_cc16_216/results.tsv');report['legacy_direct_comparisons']={'archived_method_files_checked':sum(int(v['archived_method_files_compared']) for v in z),'rows_with_differences':sum(v['all_available_archived_values_identical']!='True' for v in z)}
(cross/'final_report.json').write_text(json.dumps(report,indent=2)+'\n')
if report['core'].get('missing_archived_method_files') or report.get('legacy_direct_comparisons',{}).get('archived_method_files_checked')!=432:report['status']='REQUIRES_REVIEW'
(cross/'final_report.json').write_text(json.dumps(report,indent=2)+'\n')
lines=['# CellChatESS revision reproduction','',f"Overall status: {report['status']}.",'','## Computed outcomes','']
for x in phase:lines.append(f"- {x['stage']}: exit {x['exit_code']}, {x['elapsed_s']/3600:.2f} h")
c=report['core'];ll=report.get('legacy',{});pp=report.get('preparation_serialization',{})
lines+=['',f"Core comparison: {c.get('archived_arrays_checked',0)} same-method arrays checked against historical RDS, {c.get('different_archived_arrays','unknown')} changed, {c.get('missing_archived_method_files','unknown')} historical RDS missing. Within-run numerical acceptance failures: {c.get('fresh_acceptance_failures','unknown')}.",'',f"Legacy 1.6.1/0.1.1: {ll.get('processed',0)}/216 pairs completed. Median original/ESS inference ratio: {ll.get('pooled_median','unavailable')} (archived 70.82). The July installed CellChat build is identified by its environment; its upstream Git commit was not recorded.",'',f"Prepared inputs: {pp.get('compared_inputs',0)}/129 regenerated. {pp.get('semantic_agreement','unknown')} agreed on ordered cells and groups, normalized and signaling expression, database and LR tables. {pp.get('sha_identical','unknown')} serialized outputs also matched the frozen SHA256. The 216-to-129 ledger preserves all experimental seeds.",'','## Re-measured performance','','Times, speed ratios and peak RSS are measurements on the current server run. Runtime equality is not required.','']
for row in tsv(r/'figures/measured_claims_vs_archive.tsv') if (r/'figures/measured_claims_vs_archive.tsv').exists() else []:
 if 'CI' in row['metric']:continue
 lines.append(f"- {row['metric']}: {float(row['fresh']):.6g} fresh; {float(row['archived']):.6g} archived")
lines+=['','Full machine-readable statistics, all CIs and individual numerical differences are included in the adjacent tables. Figure1_fresh and FigureS1_fresh were drawn from the new run. The fixed group-display labels and interaction display order in panels D/E/J are reused from the original design.']
(cross/'REPORT.md').write_text('\n'.join(lines)+'\n')
from zipfile import ZipFile,ZIP_DEFLATED
selected=[]
for folder,globs in [(cross,['final_report.json','REPORT.md','phase_log.json','core/*.tsv','core/audit.json','preparation129/results.tsv','preparation129/completion.json','preparation129/all216_to_129.tsv','preparation129/source_inputs.tsv']), (r/'figures',['*.pdf','*.png','*.json','*.tsv','plot_data/*.tsv','official_statistics/plot_data/*.tsv','official_statistics/summary.json','official_statistics/all_numerical_comparisons.tsv']), (r/'legacy_cc16_216',['results.tsv','completion.json','preflight.txt']), (r/'auxiliary',['*status*.json','*summary*.json','*.tsv','core22_statistics/*.tsv','core22_statistics/*.json','ablation_CC22/all_measurements.tsv','ablation_CC22/numerical_comparisons.tsv','highest_retention_CC22/all_measurements.tsv','highest_retention_CC22/numerical_comparisons.tsv','xenium_target_trace/*summary*.json']), (r/'xenium',['fresh*.tsv','fresh*.json','stage_RSS_summary.tsv','comparison/*.tsv','comparison/verified_comparison.json','original/run_metadata.json','ess/run_metadata.json','original/stage_events.tsv','ess/stage_events.tsv']), (root/'orchestration/compare',['probe_groups.log','*.py','*.R']), (root/'workspace/results/new',['primary216/completion.json','autozyme192/completion.json','primary216/all_measurements.tsv','autozyme192/all_measurements.tsv','primary216/provenance.json','autozyme192/provenance.json'])]:
 for g in globs:selected.extend(p for p in folder.glob(g) if p.is_file())
files=sorted(set(selected))
dest=r/'CellChatESS_fresh_summary.zip'
tmp=dest.with_suffix('.zip.tmp')
with ZipFile(tmp,'w',ZIP_DEFLATED,compresslevel=6) as z:
 hashes=[]
 for f in files:
  rel=f.relative_to(root).as_posix();z.write(f,rel);hashes.append(hashlib.sha256(f.read_bytes()).hexdigest()+'  '+rel)
 z.writestr('SHA256SUMS.txt','\n'.join(hashes)+'\n')
with ZipFile(tmp) as z:assert z.testzip() is None
tmp.replace(dest)
(dest.with_suffix('.zip.sha256')).write_text(hashlib.sha256(dest.read_bytes()).hexdigest()+'  '+dest.name+'\n')
print(json.dumps({k:v for k,v in report.items() if k in ('status','preparation_serialization','legacy_direct_comparisons','core')},indent=2));raise SystemExit(int(report['status']!='PASS'))
