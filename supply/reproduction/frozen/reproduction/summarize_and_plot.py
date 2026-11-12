#!/usr/bin/env python3
"""Summarize fresh CC22 runs and draw scientific panel sources; never reads frozen plots."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=argparse.ArgumentParser();P.add_argument('--primary',type=Path,required=True);P.add_argument('--autozyme',type=Path,required=True);P.add_argument('--output',type=Path,required=True);P.add_argument('--allow-partial-smoke',action='store_true');A=P.parse_args()
if A.output.exists():P.error('Choose a new summary directory')
def load(root,expected):
 done=json.loads((root/'completion.json').read_text());prov=json.loads((root/'provenance.json').read_text())
 if done['status']!='COMPLETED':P.error('Unfinished/failed queue: '+str(root))
 if not A.allow_partial_smoke and (not done['full_run'] or done['method_processes_completed']!=expected):P.error('Full benchmark required: '+str(root))
 if prov['baseline']!='CellChat 2.2.0.9001':P.error('Wrong baseline')
 d=pd.read_csv(root/'all_measurements.tsv',sep='\t')
 if set(d.CellChat_version)!= {'2.2.0.9001'}:P.error('Mixed CellChat versions')
 if set(d.loc[d.method=='ess','ESS_loaded_version']) != {'0.2.2'}:P.error('Wrong ESS version')
 if not np.isfinite(d.inference_call_sec).all() or not (d.inference_call_sec>0).all():P.error('Invalid timings')
 return d,prov
primary,pp=load(A.primary,432);auto,ap=load(A.autozyme,576)
if pp['candidate_source_manifest_sha256']!=ap['candidate_source_manifest_sha256']:P.error('Main/AutoZyme used different candidate sources')
A.output.mkdir(parents=True);(A.output/'plot_data').mkdir();(A.output/'figures').mkdir()
def pairs(d,methods):
 idx=['experiment_id','dataset_id','scale','replicate','block']
 out=d.pivot(index=idx,columns='method',values='inference_call_sec').reset_index()
 if out[methods].isna().any().any():P.error('Unpaired methods')
 out['original_over_ess']=out.original/out.ess
 if 'autozyme' in methods:out['autozyme_over_ess']=out.autozyme/out.ess
 return out
p=pairs(primary,['original','ess']);q=pairs(auto,['original','autozyme','ess'])
p.to_csv(A.output/'plot_data/primary_pairs.tsv',sep='\t',index=False);q.to_csv(A.output/'plot_data/autozyme_pairs.tsv',sep='\t',index=False)
def cluster_summary(d,col,reps,seed):
 groups=[g[col].to_numpy() for _,g in d.groupby('dataset_id',sort=True)]
 rng=np.random.default_rng(seed);pooled=[];median_of_medians=[]
 ds=np.array([np.median(g) for g in groups]);n=len(groups)
 for _ in range(reps):
  ix=rng.integers(0,n,n);pooled.append(np.median(np.concatenate([groups[i] for i in ix])));median_of_medians.append(np.median(ds[ix]))
 return {'pooled_median':float(np.median(d[col])),'pooled_dataset_cluster_bootstrap_95CI':np.quantile(pooled,[.025,.975]).tolist(),'median_of_dataset_medians':float(np.median(ds)),'median_of_dataset_medians_bootstrap_95CI':np.quantile(median_of_medians,[.025,.975]).tolist(),'datasets':n,'records':len(d),'bootstrap_replicates':reps,'seed':seed,'generator':'numpy.random.default_rng / PCG64'}
summary={'status':'SMOKE_ONLY_NOT_MANUSCRIPT_RESULTS' if A.allow_partial_smoke else 'NEW_CC22_RESULTS','primary':cluster_summary(p,'original_over_ess',10000,20261026),'autozyme':cluster_summary(q,'autozyme_over_ess',100000,20261027)}
resolved=pd.read_csv(A.primary/'resolved_inputs.tsv',sep='\t')
z=p.merge(resolved[['experiment_id','input_sha256']],on='experiment_id',validate='one_to_one')
dedup=z.groupby(['dataset_id','input_sha256'],sort=True).agg(original=('original','median'),ess=('ess','median'),timing_records=('experiment_id','size'),original_over_ess=('original_over_ess','median')).reset_index()
dedup.to_csv(A.output/'plot_data/deduplicated_inputs.tsv',sep='\t',index=False)
summary['deduplicated_primary']=cluster_summary(dedup,'original_over_ess',10000,20261026)
summary['deduplicated_primary']['definition']='For each dataset/input SHA256, take the median of paired Original/ESS ratios; repeated full inputs are timing repeats. Never divide independent method-time medians.'
comparisons=[]
for label,root in [('primary',A.primary),('autozyme',A.autozyme)]:
 for f in sorted(root.glob('*/block*/*_vs_original.tsv')):
  d=pd.read_csv(f,sep='\t')
  # Earlier strict comparator only wrote a table after checking these fields.
  for field in ['dimensions_equal','dimnames_equal','finite_values']:
   if field not in d.columns:d[field]=True
  d['queue']=label;d['input_id']=f.parent.parent.name;d['block']=f.parent.name;d['comparison']=f.stem;comparisons.append(d)
c=pd.concat(comparisons,ignore_index=True)
c.to_csv(A.output/'all_numerical_comparisons.tsv',sep='\t',index=False)
if c[['dimensions_equal','dimnames_equal','finite_values']].isna().any().any() or not c[['dimensions_equal','dimnames_equal','finite_values']].all().all():
 (A.output/'NUMERICAL_COMPARISON_REQUIRES_DIAGNOSIS.txt').write_text('Dimension, name or nonfinite mismatch: inspect all_numerical_comparisons.tsv. Do not claim numerical agreement.\n')
 P.error('Dimensions, names or finite values differ; numerical comparisons were saved for diagnosis.')
summary['numerical']={}
for key,g in c.groupby(['queue','comparison','output']):
 summary['numerical']['/'.join(key)]={'max_absolute_difference':float(g.max_absolute_difference.max()),'different_elements_summed_over_runs':int(g.different_elements.sum()),'all_elementwise_exact':bool((g.different_elements==0).all())}
summary['autozyme']['ess_faster_blocks']=int((q.autozyme>q.ess).sum());summary['autozyme']['autozyme_faster_blocks']=int((q.autozyme<q.ess).sum());summary['autozyme']['tied_blocks']=int((q.autozyme==q.ess).sum())
summary['median_peak_RSS_GiB']=auto.groupby('method').peak_RSS_GiB.median().to_dict()
summary['candidate_source_manifest_sha256']=pp['candidate_source_manifest_sha256']
(A.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'svg.fonttype':'none','pdf.fonttype':42})
colors={'original':'#777777','autozyme':'#B15B5B','ess':'#4D79A5'}
def save(fig,name):
 for ext in ['pdf','svg','png']:fig.savefig(A.output/'figures'/f'{name}.{ext}',bbox_inches='tight',dpi=300)
 plt.close(fig)
labels=sorted(p.dataset_id.unique());fig,ax=plt.subplots(figsize=(6.6,4.2))
ax.boxplot([p.loc[p.dataset_id==s,'original_over_ess'] for s in labels],orientation='horizontal',tick_labels=labels,whis=1.5,showfliers=True)
ax.set_xlabel('Original CellChat / ESS inference-call time');ax.set_title('CellChat 2.2.0.9001 / ESS 0.2.2');save(fig,'primary_speedup_by_dataset')
med=q.groupby('dataset_id').autozyme_over_ess.median().sort_values();med.rename('median_AutoZyme_over_ESS').to_csv(A.output/'plot_data/autozyme_dataset_medians.tsv',sep='\t')
fig,ax=plt.subplots(figsize=(6.2,4.2));ax.scatter(med.values,range(len(med)),c=colors['ess']);ax.set_yticks(range(len(med)),med.index);ax.axvline(1,color='#B15B5B',ls='--');ax.set_xlabel('AutoZyme / ESS inference-call time');ax.set_title('Same CellChat 2.2.0.9001 baseline; one thread');save(fig,'autozyme_speedup_by_dataset')
fig,ax=plt.subplots(figsize=(4.6,3.4));methods=['original','autozyme','ess'];ax.boxplot([auto.loc[auto.method==m,'peak_RSS_GiB'] for m in methods],tick_labels=methods,whis=1.5);ax.set_ylabel('Independent-process peak RSS (GiB)');save(fig,'three_method_peak_RSS')
fig,ax=plt.subplots(figsize=(4.6,3.4))
for m in methods:
 d=auto[auto.method==m];ax.scatter(d.cells,d.inference_call_sec,label=m,s=12,alpha=.5,c=colors[m])
ax.set(xscale='log',yscale='log',xlabel='Cells',ylabel='Inference function call (s)');ax.legend(frameon=False);save(fig,'three_method_inference_time')
print('Saved new summaries and four panel sources. Ablation/diagnostic/RSS figure composition follows FIGURE_REBUILD_PLAN_CN.md.')
