#!/usr/bin/env python3
"""Response Figure R1: larger bold type and compact layout; measurements unchanged."""
from pathlib import Path
import argparse, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.text import Text
import pymupdf as fitz

ap=argparse.ArgumentParser()
ap.add_argument('--delivery',type=Path,default=Path(__file__).resolve().parent)
a=ap.parse_args()
B=a.delivery.resolve(); F=B/'figures'; D=F/'plot_data'; Q=F/'qa'; Q.mkdir(exist_ok=True)
PDF=B/'output/pdf';PDF.mkdir(parents=True,exist_ok=True)
S=json.loads((B/'summary_final.json').read_text())
ORIG='#B35B62'; ESS='#356F9D'; AUTO='#737C87'
plt.rcParams.update({'font.family':'DejaVu Sans','axes.spines.top':False,
    'axes.spines.right':False,'axes.linewidth':.75,'lines.linewidth':1.1,
    'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white',
    'xtick.major.width':.65,'ytick.major.width':.65,
    'xtick.major.size':3,'ytick.major.size':3})
def read(name):return pd.read_csv(D/name,sep='\t',float_precision='round_trip')
def short(s):return s.replace('3CA_data__','').replace('normal_control__','Control ')
headings=[]
def title(ax,letter,text,sub=None):
    ax.annotate(letter,xy=(0,1),xycoords='axes fraction',xytext=(-55,20),
        textcoords='offset points',size=11,weight='bold',va='bottom')
    main=ax.annotate(text,xy=(0,1),xycoords='axes fraction',xytext=(0,20),
        textcoords='offset points',size=10,weight='bold',va='bottom')
    if sub:
        subtitle=ax.annotate(sub,xy=(0,1),xycoords='axes fraction',xytext=(0,5),
            textcoords='offset points',va='bottom',size=7,color='#555')
        headings.append((letter,main,subtitle))
def save(fig,stem):
    fig.canvas.draw()
    for text in fig.findobj(match=Text):
        text.set_fontsize(text.get_fontsize()*1.14)
        text.set_fontweight('bold')
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer();bounds=fig.get_tightbbox(renderer)
    gaps={letter:(main.get_window_extent(renderer).y0-sub.get_window_extent(renderer).y1)
          /fig.dpi*72 for letter,main,sub in headings}
    assert all(gap>=4 for gap in gaps.values()),gaps
    w,h=fig.get_size_inches()
    if bounds.x0<0 or bounds.y0<0 or bounds.x1>w or bounds.y1>h:
        fig.savefig(Q/'layout_draft.png',dpi=150)
        raise RuntimeError(f'Text outside page: {bounds.bounds}; page {w,h}')
    path=PDF/f'{stem}.pdf'
    fig.savefig(path,metadata={'Title':'Figure R1. Additional computational validation of CellChatESS',
        'Subject':'Panels A-C unchanged; panel D uses complete paired Xenium memory records'})
    doc=fitz.open(path);page=doc[0]
    (F/f'{stem}.svg').write_text(page.get_svg_image(text_as_path=False))
    page.get_pixmap(dpi=600,alpha=False).save(F/f'{stem}.png')
    page.get_pixmap(dpi=150,alpha=False).save(Q/f'{stem}_preview.png')
    (Q/'layout_check.json').write_text(json.dumps({'page_mm':[297,210],
        'tight_bbox_inches':list(bounds.bounds),'all_text_weight':'bold',
        'font_scale':1.14,'data_changes':True,'changed_panels':['D'],'title_subtitle_gap_pt':gaps},indent=2)+'\n')
    doc.close();plt.close(fig)
    print(path)

p=read('A_primary_pairs_216.tsv');q=read('B_autozyme_pairs_192.tsv')
lo,hi=S['primary']['pooled_dataset_cluster_bootstrap_95CI'];x=S['xenium']

plt.rcParams.update({'font.size':8.4,'axes.titlesize':10,'axes.labelsize':8.2,'xtick.labelsize':7.6,'ytick.labelsize':7.6,'legend.fontsize':7.6})
fig=plt.figure(figsize=(297/25.4,210/25.4));gs=fig.add_gridspec(2,2,left=.13,right=.98,top=.93,bottom=.075,hspace=.42,wspace=.22)
ax=fig.add_subplot(gs[0,0]);title(ax,'A','Original CellChat / ESS speedup','CellChat 2.2.0.9001 · ESS 0.2.2')
med=p.groupby('dataset_id').original_over_ess.median().sort_values();arrays=[p[p.dataset_id.eq(s)].original_over_ess.to_numpy() for s in med.index]
box=ax.boxplot(arrays,orientation='horizontal',tick_labels=[short(s) for s in med.index],patch_artist=True,whis=1.5,showfliers=False,medianprops={'color':ORIG,'linewidth':1.2},boxprops={'facecolor':'#e6edf4','edgecolor':ESS,'linewidth':.8},whiskerprops={'color':ESS,'linewidth':.8},capprops={'color':ESS,'linewidth':.8})
for i,v in enumerate(arrays):ax.scatter(v,i+1+np.linspace(-.18,.18,len(v)),s=8,c=ESS,alpha=.36,edgecolors='none')
ax.set(xlabel='Original CellChat / ESS inference speedup (×)',ylabel='Dataset');ax.grid(axis='x',lw=.4,color='#eee');ax.set_axisbelow(True);ax.set_xlim(0,p.original_over_ess.max()*1.03);ax.set_ylim(-1.05,12.65)
ax.set_ylabel('Dataset',labelpad=-8)
ax.text(.99,.01,f'Median {S["primary"]["pooled_median"]:.2f}×; 95% dataset-cluster CI\n{lo:.2f}–{hi:.2f}×',transform=ax.transAxes,ha='right',va='bottom',size=7.5,bbox={'facecolor':'white','edgecolor':'none','pad':1.5})
ax=fig.add_subplot(gs[0,1]);title(ax,'B','AutoZyme / ESS by dataset','CellChat 2.2.0.9001 · AutoZyme 0.3.1 · ESS 0.2.2')
dm=q.groupby('dataset_id').autozyme_over_ess.median().sort_values();ax.bar(range(12),dm,color=ESS,width=.7);ax.axhline(1,color=ORIG,ls='--',lw=.9);ax.set_xticks(range(12),[short(s) for s in dm.index]);plt.setp(ax.get_xticklabels(),rotation=42,ha='right',rotation_mode='anchor',fontsize=7.1);ax.set(ylabel='Dataset median speedup (×)',xlabel='',ylim=(0,dm.max()*1.33));ax.set_xlim(-.6,11.6)
ci=S['autozyme']['median_of_dataset_medians_bootstrap_95CI'];ax.text(.02,.97,f'Median of dataset medians {S["autozyme"]["median_of_dataset_medians"]:.3f}×\n95% dataset-bootstrap CI {ci[0]:.3f}–{ci[1]:.3f}×',transform=ax.transAxes,va='top',size=7.5)
ax=fig.add_subplot(gs[1,0]);title(ax,'C','Pruning ablation: A4 vs A5','CellChat 2.2.0.9001 · ESS 0.2.2-derived variants')
hh=read('R1C_high_retention_pairs.tsv');xx=np.arange(3);ax.bar(xx-.18,hh.A4,width=.36,color=ORIG,label='A4: cached, no pruning');ax.bar(xx+.18,hh.A5,width=.36,color=ESS,label='A5: full ESS')
for i,r in hh.iterrows():ax.text(i,max(r.A4,r.A5)+.07,f'{r.A4_over_A5:.2f}×',ha='center',size=8)
ax.set_xticks(xx,['Paired run 1','Paired run 2','Paired run 3']);ax.set(xlim=(-.6,2.6),ylim=(0,hh[['A4','A5']].max().max()*1.39),ylabel='Inference-call time (s)');ax.legend(loc='upper right',frameon=False,fontsize=7.4)
ax.text(.02,.97,f'Median speedup {S["high_retention_median_speedup"]:.2f}×\nRetained work {hh.f_LR.iloc[0]*100:.3f}%',transform=ax.transAxes,va='top',size=7.8)

sub=gs[1,1].subgridspec(2,1,height_ratios=[1,1],hspace=.48)
paired_sources=[]
for i,(method,label,color) in enumerate([
    ('original','Original CellChat',ORIG),('ess','ESS',ESS)]):
    folder=EVIDENCE/method
    raw=pd.read_csv(folder/'rss_samples.tsv',sep='\t',float_precision='round_trip')
    stages=pd.read_csv(folder/'stage_events.tsv',sep='\t',float_precision='round_trip')
    run=json.loads((folder/'run_metadata.json').read_text())
    zero=run['function_call_start_epoch']*1000
    raw['elapsed_from_inference_min']=(raw.sample_start_epoch_ms-zero)/60000
    raw['RSS_GiB']=pd.to_numeric(raw.VmRSS_kB,errors='coerce')/1024**2
    stages['elapsed_from_inference_min']=(stages.timestamp_epoch_ms-zero)/60000
    ev=stages.set_index('stage')
    infer_end=float(ev.loc['after_inference','elapsed_from_inference_min'])
    agg_end=float(ev.loc['after_aggregateNet','elapsed_from_inference_min'])
    peak=float(resource_rows[i]['GNU_time_peak_RSS_GiB'])
    ax=fig.add_subplot(sub[i])
    if i==0:
        title(ax,'D','Xenium paired RSS',
              'CellChat 2.2.0.9001 · ESS 0.2.2 · single thread')
    ax.plot(raw.elapsed_from_inference_min,raw.RSS_GiB,color=color,lw=.85)
    ax.axvspan(infer_end,agg_end,color='#777777',alpha=.08,linewidth=0)
    ax.axvline(infer_end,color='#555555',ls='--',lw=.75)
    ax.axvline(agg_end,color='#555555',ls=':',lw=.9)
    ax.set(xlim=(min(-.3,float(raw.elapsed_from_inference_min.min())),
                 float(raw.elapsed_from_inference_min.max())*1.015),
           ylim=(0,280),ylabel='RSS (GiB)')
    ax.set_yticks([0,100,200])
    ax.tick_params(labelsize=7)
    ax.text(.02,.95,f'{label}: GNU peak {peak:.2f} GiB',
            transform=ax.transAxes,size=7.2,va='top',
            bbox={'facecolor':'white','edgecolor':'none','pad':.6})
    ax.set_xlabel('Elapsed from inference entry (min)',fontsize=7.1,labelpad=1)
    raw['method']=method
    raw.to_csv(D/f'R1D_{method}_paired_RSS.tsv',sep='\t',index=False)
    stages.to_csv(D/f'R1D_{method}_paired_stages.tsv',sep='\t',index=False)
    paired_sources.append({'method':method,'input_sha256':run['input']['input_sha256'],
        'inference_sec':run['inference_function_elapsed_sec'],
        'aggregate_sec':run['aggregateNet_elapsed_sec'],
        'GNU_peak_RSS_GiB':peak,
        'inference_end_min':infer_end,'aggregate_end_min':agg_end})
save(fig,'FigureR1_paired_CC22_20261102')
(Q/'paired_sources.json').write_text(json.dumps(paired_sources,indent=2)+'\n')
