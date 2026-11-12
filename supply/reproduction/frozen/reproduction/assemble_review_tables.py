#!/usr/bin/env python3
"""Read completed records only; create review tables without rerunning inference."""
from pathlib import Path
import argparse, json, shutil, subprocess, os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, pearsonr

ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
OUT=args.output.resolve();OUT.mkdir(parents=True,exist_ok=True)
T=OUT/'supplementary';D=OUT/'figures/plot_data'
T.mkdir(exist_ok=True);D.mkdir(parents=True,exist_ok=True)
N=ROOT/'results/new';S=N/'core22_statistics'
def read(p):return pd.read_csv(p,sep='\t')
def write(df,p):p.parent.mkdir(parents=True,exist_ok=True);df.to_csv(p,sep='\t',index=False)
def table(df,name,compact=None):
    write(df,T/(name+'.tsv'))
    show=df if compact is None else compact
    (T/(name+'.md')).write_text(show.to_markdown(index=False,floatfmt='.6g')+'\n')

p=read(S/'plot_data/primary_pairs.tsv');q=read(S/'plot_data/autozyme_pairs.tsv')
pm=read(N/'primary216_merged/all_measurements.tsv');am=read(N/'autozyme192/all_measurements.tsv')
f=read(S/'dataset_characteristics_216.tsv');summary=json.loads((S/'summary.json').read_text())
assert len(p)==216 and len(pm)==432 and len(q)==192 and len(am)==576
assert pm.groupby('experiment_id').method.nunique().eq(2).all()
assert am.groupby(['experiment_id','block']).method.nunique().eq(3).all()
assert am.groupby('experiment_id').block.nunique().eq(6).all()
assert pm.CellChat_version.eq('2.2.0.9001').all() and am.CellChat_version.eq('2.2.0.9001').all()
assert pm.loc[pm.method=='ess','ESS_loaded_version'].eq('0.2.2').all()
assert am.loc[am.method=='autozyme','AutoZyme_loaded_version'].eq('0.3.1').all()
assert pm.loc[pm.method=='ess','threads_used'].eq(1).all() and am.loc[am.method=='ess','threads_used'].eq(1).all()
assert np.isclose(np.median(p.original/p.ess),summary['primary']['pooled_median'])
assert np.isclose(np.median(q.autozyme/q.ess),summary['autozyme']['pooled_median'])
f['repeated_full_timing']=f.scale.astype(str).ne('all') & f.cells.lt(pd.to_numeric(f.scale,errors='coerce'))
# "all" records also belong to repeated full-input combinations. Identify by
# the same full-object hash within each dataset, rather than target labels alone.
full_hash=set(zip(f.loc[f.scale.astype(str).eq('all'),'dataset_id'],f.loc[f.scale.astype(str).eq('all'),'input_sha256']))
f['is_full_input']=pd.Series([(a,b) in full_hash for a,b in zip(f.dataset_id,f.input_sha256)])
f['records_with_same_input']=f.groupby(['dataset_id','input_sha256']).experiment_id.transform('size')
f['repeated_full_timing']=f.is_full_input & f.records_with_same_input.gt(1)
assert len(f[['dataset_id','input_sha256']].drop_duplicates())==129
assert f.positive_support_but_zero_probability.eq(0).all() and f.nonzero_probability_without_support.eq(0).all()
cols=['experiment_id','dataset_id','scale','replicate','input_sha256','cells','K','nLR','min_group_cells','median_group_cells','max_group_cells','signaling_genes','matrix_nonzero','signaling_gene_cell_elements','input_R_object_bytes','matrix_sparsity','positive_L_Radj_triples','positive_group_pair_union','f_LR','f_union','is_full_input','repeated_full_timing','records_with_same_input','original_over_ess','CellChat','CellChatESS']
table(f[cols],'S_A_dataset_characteristics_216',f.groupby('dataset_id').agg(cells_min=('cells','min'),cells_max=('cells','max'),K_min=('K','min'),K_max=('K','max'),f_LR_min=('f_LR','min'),f_LR_max=('f_LR','max')).reset_index())
ds=f.groupby('dataset_id').agg(K_min=('K','min'),K_max=('K','max'),median_speedup=('original_over_ess','median')).reset_index()
full=f[f.scale.astype(str).eq('all')].groupby('dataset_id').agg(full_cells=('cells','max'),full_K=('K','median')).reset_index()
ds=ds.merge(full,on='dataset_id',validate='one_to_one');table(ds,'S_B_dataset_summary_12')
sc=[]
for scale in ['1000','5000','10000','25000','50000','all']:
    v=p[p.scale.astype(str).eq(scale)]
    sc.append(dict(scale=scale,n=len(v),median_paired_speedup=v.original_over_ess.median(),median_saved_sec=(v.original-v.ess).median(),median_original_sec=v.original.median(),median_ess_sec=v.ess.median()))
table(pd.DataFrame(sc),'S_B2_scale_summary')
historical=read(ROOT/'reference_materials/03_results/benchmark_216/common_external_timer_216.tsv')
seeds=f[['experiment_id','dataset_id','scale','replicate','input_sha256','cells','K','is_full_input','repeated_full_timing']].merge(historical[['experiment_id','downsample_seed','seed_use']],on='experiment_id',validate='one_to_one')
assert len(seeds)==216 and seeds.downsample_seed.notna().all() and seeds.seed_use.eq(1).all()
seeds['source']='historical prepared-input seed ledger; inputs reused by SHA256, not newly downsampled'
table(seeds,'S_A2_prepared_input_seeds',seeds[['dataset_id','replicate','downsample_seed','seed_use']].drop_duplicates())
ab=read(N/'ablation_CC22/all_measurements.tsv');ac=read(N/'ablation_CC22/numerical_comparisons.tsv')
med=ab[ab.phase.eq('timing')].groupby(['dataset_id','scale','stage','cells','K','nLR'],as_index=False).agg(median_sec=('inference_call_sec','median'),min_sec=('inference_call_sec','min'),max_sec=('inference_call_sec','max'),median_peak_RSS_GiB=('peak_RSS_GiB','median'),n=('inference_call_sec','size'))
acc=ac.groupby(['experiment_id','stage','output'],as_index=False).agg(max_abs_difference=('max_absolute_difference','max'),different_elements=('different_elements','sum'))
table(med,'S_C_ablation_stage_medians');table(ab,'S_C2_ablation_all_measurements',ab[['stage','phase','replicate','scale','inference_call_sec','peak_RSS_GiB']])
table(acc,'S_C3_ablation_numerical_agreement');write(ac,T/'S_C4_ablation_all_comparisons.tsv')
ratio=[]
for scale,g in med.groupby('scale'):
    a=g.set_index('stage').median_sec
    for x,y in [('A0','A1'),('A1','A2'),('A2','A4'),('A3','A4'),('A4','A5'),('A2','A5'),('A1','A5'),('A0','A5')]:
        ratio.append(dict(scale=scale,ratio=f'{x}/{y}',ratio_of_stage_medians=a[x]/a[y]))
table(pd.DataFrame(ratio),'S_C5_ablation_ratios')
hi=ab[ab.phase.eq('high_retention_timing')].pivot(index='replicate',columns='stage',values='inference_call_sec').reset_index()
hi['A4_over_A5']=hi.A4/hi.A5
hi['f_LR']=f.loc[f.experiment_id.eq(ab.loc[ab.phase.eq('high_retention_timing'),'experiment_id'].iloc[0]),'f_LR'].iloc[0]
hi['dataset_id']='3CA_data__SmartSeq2';hi['cells']=1000;hi['K']=13
table(hi,'S_D2_earlier_rep1_high_retention_pairs')
# The handoff's nominal highest-retention row selected prepared replicate 1.
# Preserve it, and use the explicit independently recorded replicate-2 addendum.
highest=N/'highest_retention_CC22'
assert json.loads((highest/'completion.json').read_text())['status']=='COMPLETED'
ha=read(highest/'all_measurements.tsv')
hi=ha[ha.phase.eq('timing')].pivot(index='replicate',columns='stage',values='inference_call_sec').reset_index()
hi['A4_over_A5']=hi.A4/hi.A5
selected=read(highest/'selected_input.tsv').iloc[0]
assert selected.f_LR==f.f_LR.max()
hi['f_LR']=selected.f_LR;hi['dataset_id']=selected.dataset_id;hi['cells']=selected.cells;hi['K']=selected.K
hi['experiment_id']=selected.experiment_id
table(hi,'S_D_high_retention_pairs')
write(read(highest/'numerical_comparisons.tsv'),T/'S_D3_highest_retention_numerical.tsv')
inputs=am.drop_duplicates('experiment_id')[['experiment_id','dataset_id','scale','cells','K','nLR','input_sha256']]
table(inputs,'S_E_autozyme_inputs_32')
ad=q.groupby('dataset_id').agg(blocks=('block','size'),median_AutoZyme_over_ESS=('autozyme_over_ess','median'),median_Original_over_ESS=('original_over_ess','median')).reset_index()
rss=am.pivot_table(index='dataset_id',columns='method',values='peak_RSS_GiB',aggfunc='median').add_suffix('_median_RSS_GiB').reset_index()
table(ad.merge(rss,on='dataset_id'),'S_E2_autozyme_dataset_summary')
guard=[];struct=[]
for ver in ['16','22']:
    struct.append(read(N/f'structure_cc{ver}/support_pruning_matrix_results.tsv'))
    guard.append(read(N/f'input_guards_cc{ver}.tsv'))
table(pd.concat(struct,ignore_index=True),'S_F_structure_tests_52',pd.concat(struct,ignore_index=True)[['CellChat_version','structure','state','storage','probability_max_abs_diff','pvalue_different_elements','count_different_elements']])
table(pd.concat(guard,ignore_index=True),'S_F2_input_guards_10')
targets=[]
for ver in ['16','22']:
    v=read(N/f'xenium_target_trace/cc{ver}/output/target_summary.tsv');v['CellChat']='1.6.1' if ver=='16' else '2.2.0.9001';targets.append(v)
table(pd.concat(targets,ignore_index=True),'S_G_Xenium_target_summary')
x=json.loads((N/'xenium_profile/completion.json').read_text());samples=read(N/'xenium_profile/rss_samples.tsv');events=read(N/'xenium_profile/stage_events.tsv')
origin=float(events.set_index('stage').loc['R_compute_enter','epoch_ms'])
samples['elapsed_from_inference_min']=(samples.sample_start_epoch_ms-origin)/60000
samples['RSS_GiB']=samples.VmRSS_kB/1024**2
events['elapsed_from_inference_min']=(events.epoch_ms-origin)/60000
delta=np.diff(samples.sample_start_epoch_ms)/1000
x['sample_count']=len(samples);x['sampled_peak_RSS_GiB']=float(samples.RSS_GiB.max());x['max_sampling_gap_sec']=float(max(delta));x['median_sampling_interval_sec']=float(np.median(delta));x['sampling_gaps_over_2_sec']=int(sum(delta>2));x['time_zero']='Immediately before inference function call'
resource=[dict(resource=k,value=v,source='results/new/xenium_profile/completion.json' if k not in ['sample_count','sampled_peak_RSS_GiB','max_sampling_gap_sec','median_sampling_interval_sec','sampling_gaps_over_2_sec','time_zero'] else 'results/new/xenium_profile/rss_samples.tsv; stage_events.tsv') for k,v in x.items()]
table(pd.DataFrame(resource),'S_H_Xenium_resources');write(events,T/'S_H2_stage_events.tsv');write(samples,T/'S_H3_RSS_samples.tsv')
corr=[]
for col in ['matrix_nonzero','signaling_gene_cell_elements','input_R_object_bytes','positive_L_Radj_triples','f_LR']:
    rho=spearmanr(f[col],f.original_over_ess).statistic
    corr.append(dict(analysis='Spearman, records',feature=col,control='',n=len(f),rho=rho,interpretation='Descriptive; scales and timing repeats nested in datasets'))
dm=f.groupby('dataset_id')[['f_LR','original_over_ess']].median()
corr.append(dict(analysis='Spearman, dataset medians',feature='f_LR',control='',n=12,rho=spearmanr(dm.f_LR,dm.original_over_ess).statistic,interpretation='12 dataset-level medians'))
for control in ['input_R_object_bytes','matrix_nonzero','signaling_gene_cell_elements']:
    xx=rankdata(f.f_LR);yy=rankdata(f.original_over_ess);zz=rankdata(f[control]);design=np.column_stack([np.ones(len(zz)),zz])
    rx=xx-design@np.linalg.lstsq(design,xx,rcond=None)[0];ry=yy-design@np.linalg.lstsq(design,yy,rcond=None)[0]
    corr.append(dict(analysis='Partial rank, records',feature='f_LR',control=control,n=216,rho=pearsonr(rx,ry).statistic,interpretation='Pearson correlation of average-rank residuals after intercept + control rank'))
    corr.append(dict(analysis='Spearman, feature and scale',feature='f_LR',control=control,n=216,rho=spearmanr(f.f_LR,f[control]).statistic,interpretation='Association of retained fraction with the scale covariate'))
table(pd.DataFrame(corr),'S_J_feature_correlations')
write(p,D/'A_primary_pairs_216.tsv');write(q,D/'B_autozyme_pairs_192.tsv');write(med,D/'C_ablation_stage_medians.tsv');write(hi,D/'R1C_high_retention_pairs.tsv');write(samples,D/'R1D_RSS_samples.tsv');write(events,D/'R1D_stage_events.tsv')
mem=am.pivot(index=['experiment_id','dataset_id','scale','block'],columns='method',values='peak_RSS_GiB').reset_index();mem['ess_minus_autozyme_GiB']=mem.ess-mem.autozyme
write(mem,D/'GH_peak_RSS_pairs_192.tsv');write(ad,D/'R1B_dataset_medians.tsv');write(f[cols],D/'features_216.tsv')
bio=ROOT/'reference_materials/04_figures/revision_evidence_snapshot'
for name in ['J_bubble_scores.tsv','J_preserved_LR_order.tsv','J_preserved_sender_order.tsv','DE_groups.tsv','Figure1D_original_weights.tsv','Figure1E_ess_weights.tsv','biology_preservation_checks.tsv']:
    shutil.copy2(bio/'plot_data'/name,D/name)
(OUT/'figures/biology_archived').mkdir(exist_ok=True)
for name in ['Figure1D_original.pdf','Figure1E_ess.pdf']:
    shutil.copy2(bio/'panels/biology_print'/name,OUT/'figures/biology_archived'/name)
rep='3CA_data__10X__cells-1000__rep-1__engine-both__ablation-full'
of=D/'F_representative_probabilities.tsv'
if not of.exists():
    rt=json.loads((ROOT/'server_install/runtime.json').read_text())
    cmd=[rt['rscript'],'--vanilla',str(ROOT/'reproduction/export_probability_plot.R'),str(N/'primary216_merged'/rep/'block1/original/result.rds'),str(N/'primary216_merged'/rep/'block1/ess/result.rds'),str(of)]
    env={**os.environ,'R_LIBS':rt['dependency_library']};subprocess.run(cmd,env=env,check=True)
fr=read(of);summary['representative']={'input':rep,'max_abs_difference':float(fr.difference.abs().max()),'n_elements':len(fr),'nonzero_original':int(fr.original.ne(0).sum()),'K':12,'cells':1000}
summary['xenium']=x;summary['high_retention_median_speedup']=float(hi.A4_over_A5.median());summary['memory_ESS_lower_than_AutoZyme_blocks']=int(mem.ess.lt(mem.autozyme).sum())
summary['ablation_numerical']={k:{'max_abs_difference':float(g.max_absolute_difference.max()),'different_elements':int(g.different_elements.sum())} for k,g in ac.groupby('output')}
summary['feature_ranges']={c:{'min':float(f[c].min()),'max':float(f[c].max())} for c in ['cells','K','nLR','matrix_sparsity','f_LR','f_union','min_group_cells','median_group_cells','max_group_cells']}
summary['descriptive_correlations']=corr
summary['record_audit']={'primary_pairs':len(p),'primary_method_processes':len(pm),'autozyme_blocks':len(q),'autozyme_method_processes':len(am),'distinct_primary_inputs':f[['dataset_id','input_sha256']].drop_duplicates().shape[0],'repeated_full_timing_records':int(f.repeated_full_timing.sum()),'full_dataset_scale_combinations':f[f.repeated_full_timing][['dataset_id','scale']].drop_duplicates().shape[0]}
(OUT/'summary_final.json').write_text(json.dumps(summary,indent=2)+'\n')
write(pd.concat([pm.assign(queue='primary216'),am.assign(queue='autozyme192')],ignore_index=True),T/'S_K_all_main_measurements_1008.tsv')
write(read(S/'all_numerical_comparisons.tsv'),T/'S_K2_numerical_comparisons.tsv')
print('Tables and actual saved representative probability export complete:',OUT)
