#!/usr/bin/env python3
"""Read-only repository checks. No inference, no data mutation, no external calls."""
from pathlib import Path
import hashlib,json,re,ast,csv,runpy
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[1]
errors=[]
def check(ok,msg):
 if not ok:errors.append(msg)
files=[p for p in R.rglob('*') if p.is_file() and '.git' not in p.parts and '__pycache__' not in p.parts]
check(all(p.stat().st_size<=50_000_000 for p in files),'File over 50MB')
check(not any(p.suffix.lower()=='.rds' for p in files),'RDS binary bundled')
secret=re.compile(r'(?:ghp_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}|sk-(?:proj-)?[A-Za-z0-9_-]{25,}|-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----)')
private=re.compile(r'/(?:home|data|Users|mnt|tmp|root)/|\b(?:dzf|alan)\b')
# Match literal private roots in data/docs, excluding this regex-based checker itself.
for p in files:
 if p==Path(__file__).resolve() or p.suffix.lower() in ['.pdf','.png']:continue
 t=p.read_text(errors='replace')
 check(not secret.search(t),'Credential-pattern match: '+str(p.relative_to(R)))
 check(not private.search(t),'Private-path/username match: '+str(p.relative_to(R)))
 if p.suffix=='.py':
  try:ast.parse(t)
  except SyntaxError:errors.append('Python parse error '+str(p.relative_to(R)))
a=pd.read_csv(R/'results/benchmark_216/common_external_timer_216.tsv',sep='\t');b=pd.read_csv(R/'results/autozyme_192/paired_inference_192.tsv',sep='\t');f=pd.read_csv(R/'results/dataset_characteristics/dataset_characteristics_216.tsv',sep='\t');p1=pd.read_csv(R/'results/benchmark_216/P1_accuracy_216.tsv',sep='\t');p2=pd.read_csv(R/'results/benchmark_216/representative/P2_representative_accuracy.tsv',sep='\t')
check(len(a)==216 and a.dataset_id.nunique()==12,'Main216 counts')
check(np.allclose(a.paired_speedup,a.original_inference_call_sec/a.ess_inference_call_sec),'Main paired ratios')
check(round(a.paired_speedup.median(),2)==70.82,'Main70.82')
check(len(b)==192 and round(b.autozyme_over_ess.median(),3)==1.669,'Compatible counts/median')
check((b.autozyme_over_ess>1).sum()==183,'ESS183blocks')
check(len(f)==216 and f.repeated_full_input_timing.sum()==99,'Repeated full99')
check(p1.pvalue_unequal_elements.sum()==0 and p1.aggregate_count_exact.all(),'Main accuracy')
check(p2.loc[p2.comparison=='new_ESS_vs_saved_sparse_stream','max_abs_prob_diff'].iloc[0]==0,'Representative exact saved output')
support=pd.read_csv(R/'results/dataset_characteristics/F5_support_set_counts_216.tsv',sep='\t');check(len(support)==216 and support.sets_equal.all(),'Support sets')
timing=pd.read_csv(R/'results/TIMING_LOG_INDEX.tsv',sep='\t',keep_default_na=False)
timefiles={str(p.relative_to(R)) for p in files if p.name.endswith('.time.txt') or p.name=='GNU_time.txt'}
check(len(timing)==584 and timing.log_path.is_unique and set(timing.log_path)==timefiles,'Timing index coverage')
for z in timing.itertuples():
 check((R/z.record_path).is_file(),'Timing-index record missing '+z.record_path)
 check(hashlib.sha256((R/z.log_path).read_bytes()).hexdigest()==z.log_sha256,'Timing log changed '+z.log_path)
 check(all(str(getattr(z,k)).strip() for k in ['experiment_id','input_id','method','CellChat','CellChatESS','AutoZyme']),'Timing-index fields missing '+z.log_path)
check((R/'supplementary/S_I_versions.tsv').read_bytes()==(R/'supplementary/tables/S_I_versions.tsv').read_bytes(),'Table R3 alias differs from S_I')
panel=pd.read_csv(R/'figures/panel_environment_map.tsv',sep='\t').set_index('file').environment
check(panel['figures/plot_data/FigureR1_supplementary_final_A_box_statistics.tsv']=='main: CellChat 1.6.1 + CellChatESS 0.1.1','R1A environment')
check('not a released' in panel['figures/plot_data/R1C_high_retention_pairs.tsv'],'A4 intermediate build label')
derived=pd.read_csv(R/'results/xenium/bootstrap_comparison_counts_216_derived.tsv',sep='\t')
check(len(derived)==216 and int(derived.bootstrap_comparisons_derived.sum())==413671500,'216 derived bootstrap total')
check((derived.bootstrap_comparisons_derived==derived.positive_support_triples*derived.nboot).all(),'216 derived bootstrap formula')
negroot=R/'results/pruning_tests/negative_ligand_evidence'
negative=pd.read_csv(negroot/'pruning_boundary_checks_original.tsv',sep='\t',float_precision='round_trip')
baseline=negative[(negative['case']=='negative_ligand')&(negative.engine=='original')]
check(len(baseline)==1 and np.isclose(baseline.prob_min.iloc[0],-0.000279975558514351,rtol=1e-14,atol=0),'Negative ligand measured minimum')
negative_provenance=pd.read_csv(negroot/'SOURCE_PROVENANCE.tsv',sep='\t')
source=negative_provenance[negative_provenance.copy_kind=='byte-identical original result'].iloc[0]
check(hashlib.sha256((R/source.copied_file).read_bytes()).hexdigest()==source.source_sha256,'Negative-ligand original copy hash')
for code,output in [
 ('results/ablation/recompute_Q3_from_indices.py','results/ablation/Q3_recomputed_from_exported_indices.tsv'),
 ('results/xenium/recompute_Q7_comparison_counts.py','results/xenium/Q7_recomputed_comparison_summary.tsv'),
]:
 try:
  # run_path defines the pure calculate function; __main__ output writers do not run.
  recalculated=runpy.run_path(str(R/code))['calculate']()
  with (R/output).open() as stream:saved=list(csv.DictReader(stream,delimiter='\t'))
  check(saved==[{key:str(value) for key,value in row.items()} for row in recalculated],'Q3/Q7 saved summary differs from fresh calculation '+output)
 except Exception as exc:errors.append('Q3/Q7 source verification failed '+code+': '+repr(exc))
provenance=pd.read_csv(R/'results/Q3_Q7_SOURCE_PROVENANCE.tsv',sep='\t')
for z in provenance.itertuples():
 check(hashlib.sha256((R/z.repository_path).read_bytes()).hexdigest()==z.repository_sha256,'Q3/Q7 evidence copy changed '+z.repository_path)
 if z.copy_kind=='byte-identical copy':check(z.source_sha256==z.repository_sha256,'Q3/Q7 raw copy not byte identical '+z.repository_path)
qchecks=pd.read_csv(R/'results/number_consistency_check_Q3_Q7.tsv',sep='\t')
check(len(qchecks)==4 and qchecks.status.eq('matches').all() and qchecks.actual.eq(qchecks.expected).all(),'Q3/Q7 numeric consistency')
manifest=R/'MANIFEST.tsv'
if manifest.exists():
 m=pd.read_csv(manifest,sep='\t');actual={str(p.relative_to(R)) for p in files if p!=manifest};check(set(m.path)==actual,'Manifest file set differs')
 for z in m.itertuples():
  p=R/z.path
  check(p.exists() and p.stat().st_size==z.size_bytes and hashlib.sha256(p.read_bytes()).hexdigest()==z.sha256,'Manifest mismatch '+z.path)
known=pd.read_csv(R/'results/number_consistency_check.tsv',sep='\t');check(len(known[known.status=='MISMATCH'])==0,'Unexpected number mismatches')
ab=known[known.metric=='A1/A5 1k'].iloc[0];check(round(float(ab.observed),2)==2.59 and ab.status=='matches','A1/A5 1k approved rounding')
print(json.dumps({'files_checked':len(files),'errors':errors,'declared_number_mismatches':known[known.status=='MISMATCH'].to_dict('records'),'pending_evidence':[],'Q3_Q7_checks':'recomputed from exported indices and independent saved support sources'},indent=2))
raise SystemExit(bool(errors))
