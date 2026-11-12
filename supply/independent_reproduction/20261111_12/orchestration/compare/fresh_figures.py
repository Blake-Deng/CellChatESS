"""Replot A–J and S1 from the fresh campaign. Only the a priori label/selection order is read from archived plotting inputs."""
import argparse,csv,json,sys,subprocess,shutil
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from scipy.stats import rankdata,spearmanr
def load(path):return pd.read_csv(path,sep='\t',float_precision='round_trip')
def write(tab,path):tab.to_csv(path,sep='\t',index=False,float_format='%.17g')
def main():
 ap=argparse.ArgumentParser(__doc__);ap.add_argument('--root',required=True,type=Path);a=ap.parse_args();r=a.root.resolve();w=r/'workspace';sup=r/'repository/supply';pr=w/'results/new/primary216';az=w/'results/new/autozyme192';aux=r/'results/auxiliary';xen=r/'results/xenium';fig=r/'results/figures';fig.mkdir(exist_ok=False,parents=True);dt=fig/'plot_data';dt.mkdir()
 rs=json.loads((w/'server_install/runtime.json').read_text())['rscript'];frozen=w/'reproduction';old=sup/'legacy/main_figure_20261101_requires_Xenium_I_update';sty=old/'plot_data'
 subprocess.run([sys.executable,str(frozen/'summarize_and_plot.py'),'--primary',str(pr),'--autozyme',str(az),'--output',str(fig/'official_statistics')],check=True)
 off=fig/'official_statistics';stats=json.loads((off/'summary.json').read_text());p=load(off/'plot_data/primary_pairs.tsv');q=load(off/'plot_data/autozyme_pairs.tsv');mem=load(az/'all_measurements.tsv').pivot(index=['experiment_id','dataset_id','scale','replicate','block'],columns='method',values='peak_RSS_GiB').reset_index();mem['saved_MiB']=(mem.autozyme-mem.ess)*1024
 assert (len(p),len(q),len(mem))==(216,192,192)
 write(p,dt/'A_primary.tsv');write(q,dt/'B_AutoZyme.tsv');write(mem,dt/'G_H_memory.tsv')
 feature=load(aux/'core22_statistics/dataset_characteristics_216.tsv').drop(columns=['original_over_ess','original','ess'],errors='ignore').merge(p[['experiment_id','original_over_ess','original','ess']],on='experiment_id',validate='one_to_one');feature['full_cells']=feature.groupby('dataset_id').cells.transform('max');feature['is_full_input']=feature.cells==feature.full_cells;feature['records_with_same_input']=feature.groupby('input_sha256').input_sha256.transform('size');feature['repeated_full_timing']=feature.is_full_input & (feature.records_with_same_input>1)
 if not 'ESS_loaded_version' in feature:feature['CellChatESS']='0.2.2'
 tmpl=load(sup/'data/Supplementary_Data_S7b.tsv');ordered=feature.set_index('experiment_id').loc[tmpl.experiment_id].reset_index();write(ordered[list(tmpl.columns)],dt/'Supplementary_Data_S7b_fresh.tsv')
 f=feature; ranks=lambda x:rankdata(x,method='average');rho=lambda x,y:float(spearmanr(x,y).statistic)
 residual=lambda x,c:ranks(x)-np.column_stack([np.ones(len(x)),ranks(c)])@np.linalg.lstsq(np.column_stack([np.ones(len(x)),ranks(c)]),ranks(x),rcond=None)[0]
 gd=f.groupby('dataset_id')[['f_LR','original_over_ess']].median()
 stats['feature_associations']={'rho_speedup_nonzero':rho(f.original_over_ess,f.matrix_nonzero),'rho_speedup_fLR':rho(f.original_over_ess,f.f_LR),'rho_dataset_medians':rho(gd.original_over_ess,gd.f_LR),'partial_fLR_rank_given_nonzero':float(np.corrcoef(residual(f.original_over_ess,f.matrix_nonzero),residual(f.f_LR,f.matrix_nonzero))[0,1]),'fLR_range':[float(f.f_LR.min()),float(f.f_LR.max())],'f_union_range':[float(f.f_union.min()),float(f.f_union.max())]}
 f['saved_sec']=f.original-f.ess;write(f.groupby('scale').agg(paired_median_speedup=('original_over_ess','median'),median_saved_sec=('saved_sec','median')).reset_index(),dt/'scale_effect.tsv')
 ab=load(aux/'ablation_CC22/all_measurements.tsv').query("phase=='timing'");ab=ab.groupby(['scale','dataset_id','stage']).agg(median_sec=('inference_call_sec','median'),peak_rss_GiB=('peak_RSS_GiB','median'),repeats=('stage','size')).reset_index();assert len(ab)==18 and ab.repeats.eq(3).all();write(ab,dt/'C_ablation.tsv');abp=ab.pivot(index='scale',columns='stage',values='median_sec');ar=pd.DataFrame({x+'_over_'+y:abp[x]/abp[y] for x,y in [('A0','A1'),('A1','A2'),('A4','A5'),('A1','A5')]});write(ar.reset_index(),dt/'C_ablation_ratios.tsv')
 hi=load(aux/'highest_retention_CC22/all_measurements.tsv').query("phase=='timing'").pivot(index='replicate',columns='stage',values='inference_call_sec').reset_index();hi['A4_over_A5']=hi.A4/hi.A5;assert len(hi)==3;write(hi,dt/'S1C.tsv')
 xr=load(xen/'fresh_R2.tsv').set_index('method');write(xr.reset_index(),dt/'I_xenium.tsv')
 F_id='3CA_data__10X__cells-1000__rep-1__engine-both__ablation-full'
 N_id='3CA_data__PRJNA289244__cells-1000__rep-1__engine-both__ablation-full'
 fb=pr/F_id/'block1';nb=pr/N_id/'block1'
 subprocess.run([rs,'--vanilla',str(frozen/'export_probability_plot.R'),str(fb/'original/result.rds'),str(fb/'ess/result.rds'),str(dt/'F_scatter.tsv')],check=True)
 subprocess.run([rs,'--vanilla',str(old/'native_code/figure1/export_DEJ.R'),str(nb),str(dt)],check=True)
 groups=load(sty/'DEJ_group_display_mapping.tsv');prepared=load(pr/'resolved_inputs.tsv').set_index('experiment_id').loc[N_id]
 subprocess.run([rs,'--vanilla','-e','suppressPackageStartupMessages(library(CellChat)); a<-commandArgs(TRUE);o<-readRDS(a[1]);write.table(as.data.frame(table(o@idents)),a[2],row.names=FALSE,sep="\\t")',str(prepared.input_file),str(dt/'raw_group_counts.tsv')],check=True)
 fresh_groups=load(dt/'raw_group_counts.tsv');fresh_groups.columns=['cluster_id','cells'];gg=fresh_groups.merge(groups,on=['cluster_id','cells'],validate='one_to_one');assert len(gg)==9 and gg.cells.sum()==1000
 gg=gg.set_index('group').loc[load(sty/'DE_groups.tsv').group].reset_index();write(gg,dt/'D_E_display_groups.tsv')
 lrs=load(sty/'J_preserved_LR_order.tsv');senders=load(sty/'J_preserved_sender_order.tsv');pairs=load(sty/'J_bubble_scores.tsv')[['source','target','pair','pair_rank']].drop_duplicates();assert len(pairs)==24
 def labels(x):return x.map(gg.set_index('cluster_id').group)
 all_lr=load(dt/'DEJ_nonzero_probabilities.tsv');all_lr['source']=labels(all_lr.source_id);all_lr['target']=labels(all_lr.target_id)
 jj=all_lr[all_lr.source.isin(senders.sender)&all_lr.interaction_name.isin(lrs.interaction_name)&(all_lr.probability>0)&(all_lr.pvalue<.05)].merge(lrs,on='interaction_name',validate='many_to_one').merge(pairs,on=['source','target'],how='left',validate='many_to_one')
 if jj.pair_rank.isna().any():raise ValueError('new significant group pair outside original display map')
 jj['score']=-1/np.log(jj.probability);write(jj,dt/'J_bubbles.tsv')
 nw=load(dt/'DE_network_weights.tsv'); fr=load(dt/'F_scatter.tsv');stats['panel_F_max_absolute_difference']=float(fr.difference.abs().max());stats['panel_J_points_per_method']=jj.groupby('engine').size().to_dict();stats['xenium']=xr.to_dict(orient='index');stats['peak_memory_blocks_ESS_lower']=int((mem.saved_MiB>0).sum())
 # Fresh statistical estimates and references are labeled separately.
 hist=json.loads((sup/'results/core22_statistics/summary.json').read_text());claims=[]
 def claim(name,new,old):claims.append(dict(metric=name,fresh=float(new),archived=float(old),fresh_minus_archive=float(new)-float(old)))
 for m in ['primary','autozyme','deduplicated_primary']:
  for k in ['pooled_median','median_of_dataset_medians']:claim(m+' / '+k,stats[m][k],hist[m][k])
  for k in ['pooled_dataset_cluster_bootstrap_95CI','median_of_dataset_medians_bootstrap_95CI']:
   for i in range(2):claim(m+' / '+k+' / '+str(i),stats[m][k][i],hist[m][k][i])
 claim('ESS faster vs AutoZyme, blocks',stats['autozyme']['ess_faster_blocks'],hist['autozyme']['ess_faster_blocks'])
 for k,v in stats['median_peak_RSS_GiB'].items():claim(k+' median peak RSS, GiB',v,hist['median_peak_RSS_GiB'][k])
 for k,oldval in [('inference_sec',13.4),('wall_sec',8.03)]:claim('Xenium speedup, '+k,xr.loc['original',k]/xr.loc['ess',k],oldval)
 claim('Xenium peak RSS reduction percent',100*(1-xr.loc['ess','peak_RSS_GiB']/xr.loc['original','peak_RSS_GiB']),20.58);claim('AutoZyme/ESS memory yield count',stats['peak_memory_blocks_ESS_lower'],192);claim('Highest retention A4/A5 median',hi.A4_over_A5.median(),1.83);write(pd.DataFrame(claims),fig/'measured_claims_vs_archive.tsv')
 stats['figure_basis']='All measurements, networks, scatter, bubble probabilities and RSS traces are from this run. Only archived display labels/order define the D/E/J visual mapping.';(fig/'statistics.json').write_text(json.dumps(stats,indent=2,default=str)+'\n')
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.titlesize':8,'axes.titleweight':'bold'})
 blue,red,gray='#165B90','#AB4054','#777777';sz=['1000','5000','10000','25000','50000','all'];strata=['small_1k','medium_10k_or_full','large_50k_or_full']
 fi=plt.figure(figsize=(11,13),layout='constrained');gr=fi.add_gridspec(4,3,height_ratios=[1,1,1,1.28]);axes=[fi.add_subplot(gr[i//3,i%3]) for i in range(9)]
 def title(ax,st):ax.set_title(st,loc='left',pad=5)
 def dist(ax,data,colmap,scales,labels):
  for j,(c,color) in enumerate(colmap):
   for i,sz0 in enumerate(scales):
    v=data[data.scale.astype(str)==sz0][c].to_numpy();xx=i+(j-(len(colmap)-1)/2)*.21;assert len(v)>0
    ax.scatter(xx+np.linspace(-.045,.045,len(v)),v,s=4,alpha=.25,color=color);t0,m,t1=np.quantile(v,[.25,.5,.75]);ax.errorbar(xx,m,yerr=[[m-t0],[t1-m]],fmt='o',ms=3,color=color,capsize=2,label=c if i==0 else None)
  ax.set_xticks(range(len(scales)),labels,rotation=25);ax.legend(fontsize=6,frameon=False)
 ax=axes[0];title(ax,'A   Original CellChat vs ESS');dist(ax,p,[('original',gray),('ess',blue)],sz,['1k','5k','10k','25k','50k','all']);ax.set_yscale('log');ax.set_ylabel('Inference call (s)');ax.set_xlabel(f"Median speedup {stats['primary']['pooled_median']:.2f}×")
 ax=axes[1];title(ax,'B   AutoZyme vs ESS');ax.scatter(q.autozyme,q.ess,s=7,color=blue);lo=min(q.ess.min(),q.autozyme.min());up=max(q.ess.max(),q.autozyme.max());ax.plot([lo,up],[lo,up],':k');ax.set(xscale='log',yscale='log',xlabel='AutoZyme inference (s)',ylabel='ESS inference (s)');ax.text(.03,.96,f"{stats['autozyme']['pooled_median']:.3f}×, {stats['autozyme']['ess_faster_blocks']}/192 faster",transform=ax.transAxes,va='top')
 ax=axes[2];title(ax,'C   Six-stage ablation')
 for k,aa in ab.groupby('scale'):ax.plot(aa.stage,aa.median_sec,'-o',ms=3,label=str(k))
 ax.set_yscale('log');ax.set_ylabel('Median inference time (s)');ax.legend(fontsize=6)
 angle=np.linspace(0,2*np.pi,9,endpoint=False)+np.pi/2;xy=dict(zip(gg.cluster_id,zip(np.cos(angle),np.sin(angle))));colors=gg.set_index('cluster_id').colour
 for idx,who in [(3,'Original CellChat'),(4,'CellChatESS')]:
  ax=axes[idx];title(ax,('D   ' if idx==3 else 'E   ')+who+' | oligodendroglioma')
  for z in nw[nw.engine==who].itertuples():
   if z.weight<=0:continue
   s,t=xy[z.source_id],xy[z.target_id];same=s==t
   if same:t=(s[0]+.055,s[1]+.06)
   ax.add_patch(FancyArrowPatch(s,t,arrowstyle='-|>',mutation_scale=4,linewidth=.13+2.0*z.weight/nw.weight.max(),alpha=.45,color=colors[z.source_id],shrinkA=4,shrinkB=4,connectionstyle='arc3,rad='+('1.7' if same else '.12')))
  for z in gg.itertuples():
   x,y=xy[z.cluster_id];ax.scatter([x],[y],s=12+.6*z.cells,c=z.colour,zorder=3);ax.text(x*1.22,y*1.18,z.group.replace('Malignant-','M'),ha='center',va='center',fontsize=6)
  ax.set_aspect('equal');ax.set(xlim=(-1.55,1.55),ylim=(-1.42,1.42));ax.axis('off')
 ax=axes[5];title(ax,'F   Numerical agreement');ax.scatter(fr.original,fr.ESS,s=2,alpha=.55,rasterized=True);up=fr.original.max();ax.plot([0,up],[0,up],':k');ax.set(xlabel='Original probability',ylabel='ESS probability');ax.text(.05,.92,f'Max |Δ| = {fr.difference.abs().max():.2e}',transform=ax.transAxes)
 ax=axes[6];title(ax,'G   Peak process RSS');dist(ax,mem,[('original',gray),('autozyme',red),('ess',blue)],strata,['1k','10k/full','50k/full']);ax.set_ylabel('Peak RSS (GiB)')
 ax=axes[7];title(ax,'H   Peak RSS saved by ESS');blocks=[mem[mem.scale.eq(s)].saved_MiB for s in strata];ax.boxplot(blocks,orientation='horizontal',tick_labels=['1k','10k/full','50k/full'],whis=(0,100));ax.axvline(0,ls=':',color=gray);ax.set_xlabel('AutoZyme minus ESS (MiB)')
 ax=axes[8];title(ax,'I   Non-spatial Xenium');ax.axis('off');tab=ax.table(cellText=[[l,f'{xr.loc["original",k]:,.2f}',f'{xr.loc["ess",k]:,.2f}'] for k,l in [('inference_sec','Inference (s)'),('wall_sec','Whole process (s)'),('peak_RSS_GiB','Peak RSS (GiB)')]],colLabels=['Measure','Original','ESS'],bbox=[0,.35,1,.52]);tab.auto_set_font_size(False);tab.set_fontsize(6);ax.text(0,.22,f"Inference {xr.loc['original','inference_sec']/xr.loc['ess','inference_sec']:.2f}×; whole process {xr.loc['original','wall_sec']/xr.loc['ess','wall_sec']:.2f}×",fontsize=7)
 sub=gr[3,:].subgridspec(1,2);vlo,vhi=jj.score.min(),jj.score.max()
 for k,who in enumerate(['Original CellChat','CellChatESS']):
  ax=fi.add_subplot(sub[k]);title(ax,('J   ' if k==0 else '')+who);zz=jj[jj.engine.eq(who)];sc=ax.scatter(zz.pair_rank,zz.lr_rank,c=zz.score,s=4+28*(zz.score-vlo)/(vhi-vlo),cmap='Blues',vmin=vlo,vmax=vhi);pair=pairs.sort_values('pair_rank');ax.set_xticks(pair.pair_rank,pair.pair,rotation=72,ha='right',fontsize=5.5);ax.set_yticks(lrs.lr_rank,lrs.interaction_name if not k else ['']*len(lrs));ax.set_ylim(16.8,.3);ax.tick_params(axis='y',labelsize=6);ax.grid(alpha=.16)
  if k:
   cb=fi.colorbar(sc,ax=ax,shrink=.55);raw=np.array([.001,.01,.05,.1,.3]);valid=raw[(-1/np.log(raw)>=vlo)&(-1/np.log(raw)<=vhi)];cb.set_ticks(-1/np.log(valid),labels=[str(x) for x in valid]);cb.set_label('Communication probability')
 fi.savefig(fig/'Figure1_fresh.pdf');fi.savefig(fig/'Figure1_fresh.png',dpi=180);plt.close(fi)
 f2=plt.figure(figsize=(11,7),layout='constrained');g=f2.add_gridspec(2,2);ax=f2.add_subplot(g[0,0]);title(ax,'A   Original/ESS speedup by dataset');lab=p.groupby('dataset_id').original_over_ess.median().sort_values().index;sh=lambda s:s.replace('3CA_data__','').replace('normal_control__','');ax.boxplot([p[p.dataset_id==x].original_over_ess for x in lab],orientation='horizontal',tick_labels=[sh(x) for x in lab]);ax.set_xlabel('Inference speedup (×)')
 ax=f2.add_subplot(g[0,1]);title(ax,'B   AutoZyme/ESS speedup by dataset');med=q.groupby('dataset_id').autozyme_over_ess.median().sort_values();ax.barh([sh(x) for x in med.index],med,color=blue);ax.axvline(1,color=gray,ls=':');ax.set_xlabel('Median inference speedup (×)')
 ax=f2.add_subplot(g[1,0]);title(ax,'C   Pruning at highest retained fraction');v=np.arange(3);ax.bar(v-.2,hi.A4,.4,label='A4',color=red);ax.bar(v+.2,hi.A5,.4,label='A5',color=blue);ax.set_xticks(v,[f'Pair {i}' for i in hi.replicate]);ax.set_ylabel('Inference time (s)');ax.legend();ax.set_xlabel(f"Median speedup {hi.A4_over_A5.median():.2f}×")
 sub=g[1,1].subgridspec(2,1)
 for i,m in enumerate(['original','ess']):
  ax=f2.add_subplot(sub[i]);rss=load(xen/m/'rss_samples.tsv');event=load(xen/m/'stage_events.tsv').set_index('stage');meta=json.loads((xen/m/'run_metadata.json').read_text());t0=meta['function_call_start_epoch']*1000;tt=(rss.sample_start_epoch_ms-t0)/60000;breaks=tt.diff()>3/60;yy=rss.VmRSS_kB/1024**2;ax.plot(tt,yy.mask(breaks),color=[gray,blue][i],lw=.75);ax.set_ylabel('RSS (GiB)');ax.set_title(('D   ' if i==0 else '')+m+'   |   Non-spatial Xenium',loc='left');t1=(event.at['after_inference','timestamp_epoch_ms']-t0)/60000;t2=(event.at['after_aggregateNet','timestamp_epoch_ms']-t0)/60000;ax.axvspan(t1,t2,color='#aaaaaa',alpha=.2);ax.axvline(t1,color=gray,ls='--',lw=.8);ax.axvline(t2,color=gray,ls=':',lw=.8);ax.set_xlabel('Minutes from inference entry')
 f2.savefig(fig/'FigureS1_fresh.pdf');f2.savefig(fig/'FigureS1_fresh.png',dpi=180);plt.close(f2)
 print('Fresh figures, source values and numerical/performance comparisons:',fig)
if __name__=='__main__':main()
