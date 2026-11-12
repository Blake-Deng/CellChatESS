#!/usr/bin/env python3
"""Figure 1 with D/E/J recomputed from saved CellChat22/ESS022 results."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec, ticker, colors
from matplotlib.lines import Line2D
from matplotlib.text import Text
import pymupdf as fitz
from network_panels import draw_network

ap=argparse.ArgumentParser();ap.add_argument('--delivery',type=Path,required=True);a=ap.parse_args()
B=a.delivery.resolve();F=B/'figures';D=F/'plot_data';Q=F/'qa';Q.mkdir(exist_ok=True)
PDF=B/'output/pdf';PDF.mkdir(parents=True,exist_ok=True)
S=json.loads((B/'summary_final.json').read_text())
ORIG='#B35B62';ESS='#356F9D';AUTO='#737C87';PALETTE=[ESS,ORIG,'#5D8B8B']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.5,'axes.titlesize':8.1,'axes.labelsize':7.2,'xtick.labelsize':6.7,'ytick.labelsize':6.7,'legend.fontsize':6.2,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.75,'lines.linewidth':1.1,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','xtick.major.width':.65,'ytick.major.width':.65,'xtick.major.size':2.5,'ytick.major.size':2.5})
def read(n):return pd.read_csv(D/n,sep='\t')
def short(s):return s.replace('3CA_data__','').replace('normal_control__','Control ')
def title(ax,letter,text,sub=None):
    ax.text(-.18,1.10,letter,transform=ax.transAxes,size=11,weight='bold',va='bottom')
    ax.set_title(text,loc='left',weight='bold',pad=15)
    if sub:ax.text(0,1.025,sub,transform=ax.transAxes,va='bottom',size=6.2,color='#555')
def log(ax,which='y'):
    getattr(ax,'set_'+which+'scale')('log');getattr(ax,which+'axis').set_major_formatter(ticker.FuncFormatter(lambda x,p:f'{x:g}'))
    getattr(ax,which+'axis').set_minor_locator(ticker.NullLocator())
def dist(ax,df,methods,scale,order,labels):
    for j,s in enumerate(order):
        for i,(method,col,color) in enumerate(methods):
            v=df[df[scale].astype(str).eq(s)][col].to_numpy();x=j+(i-(len(methods)-1)/2)*.23
            ax.scatter(x+np.linspace(-.05,.05,len(v)),v,s=3,alpha=.22,color=color,edgecolors='none')
            q,m,z=np.quantile(v,[.25,.5,.75]);ax.errorbar(x,m,yerr=[[m-q],[z-m]],fmt='o',ms=3,capsize=2,color=color,lw=1)
    ax.set_xticks(range(len(order)),labels);ax.set_xlim(-.5,len(order)-.5)
def save(fig,stem,pdfs=None):
    # Resolve every tick and table label before applying the common text style.
    fig.canvas.draw()
    for artist in fig.findobj(match=Text):
        artist.set_fontsize(artist.get_fontsize()*1.14)
        artist.set_fontweight('bold')
    for axis in fig.axes:
        for table in axis.tables:
            for cell in table.get_celld().values():
                text=cell.get_text()
                text.set_fontsize(text.get_fontsize()*1.14)
                text.set_fontweight('bold')
    fig.canvas.draw();ren=fig.canvas.get_renderer();bb=fig.get_tightbbox(ren);w,h=fig.get_size_inches()
    if bb.x0<-.015 or bb.y0<-.015 or bb.x1>w+.015 or bb.y1>h+.015:
        fig.savefig(Q/'layout_draft.png',dpi=150)
        raise RuntimeError(f'Page clipping {stem}: {bb.bounds}; page {w,h}')
    layout=Q/f'{stem}_layout.pdf';fig.savefig(layout)
    doc=fitz.open(layout);page=doc[0]
    if pdfs:
        for bounds,src in pdfs:
            x,y,w,h=bounds;rect=fitz.Rect(x*page.rect.width,(1-y-h)*page.rect.height,(x+w)*page.rect.width,(1-y)*page.rect.height)
            old=fitz.open(src);blocks=[fitz.Rect(b['bbox']) for b in old[0].get_text('dict')['blocks'] if b['type']==0]
            clip=fitz.Rect(blocks[0])
            for b in blocks[1:]:clip|=b
            clip=fitz.Rect(clip.x0-12,clip.y0-12,clip.x1+12,clip.y1+12)&old[0].rect
            page.show_pdf_page(rect,old,0,clip=clip,keep_proportion=True);old.close()
    doc.set_metadata({'title':stem,'subject':'Verified CellChat 2.2.0.9001 / ESS 0.2.2 results; biological D/E/J from saved final-candidate outputs','creator':'Matplotlib and PyMuPDF'})
    doc.save(PDF/f'{stem}.pdf',garbage=4,deflate=True)
    (F/f'{stem}.svg').write_text(page.get_svg_image(text_as_path=False))
    page.get_pixmap(dpi=600,alpha=False).save(F/f'{stem}.png')
    page.get_pixmap(dpi=150,alpha=False).save(Q/f'{stem}_preview.png')
    (Q/f'{stem}_page_check.json').write_text(json.dumps({'page_mm':[page.rect.width/72*25.4,page.rect.height/72*25.4],'page_count':len(doc),'raster_dpi':600,'tight_bbox_inches':list(bb.bounds),'text_contained':True,'vector_biological_panels':True},indent=2))
    doc.close();plt.close(fig);layout.unlink()

p=read('A_primary_pairs_216.tsv');q=read('B_autozyme_pairs_192.tsv');ab=read('C_ablation_stage_medians.tsv');mem=read('GH_peak_RSS_pairs_192.tsv');fr=read('F_representative_probabilities.tsv')
scales=['1000','5000','10000','25000','50000','all'];labs=['1k','5k','10k','25k','50k','all'];strata=['small_1k','medium_10k_or_full','large_50k_or_full'];sl=['1k','10k\nor full','50k\nor full']
fig=plt.figure(figsize=(210/25.4,270/25.4));gs=fig.add_gridspec(4,3,left=.075,right=.98,top=.95,bottom=.17,hspace=.52,wspace=.45,height_ratios=[1,1.22,1.08,1.4])
ax=fig.add_subplot(gs[0,0]);title(ax,'A','Original CellChat vs ESS','12 datasets · 216 paired runs')
dist(ax,p,[('Original','original',ORIG),('ESS','ess',ESS)],'scale',scales,labs);log(ax)
ax.set(xlabel='Target cells',ylabel='Inference-call time (s; log)');ax.set_ylim(p[['ess','original']].min().min()*.55,p.original.max()*6)
lo,hi=S['primary']['pooled_dataset_cluster_bootstrap_95CI']
ax.text(.03,.97,f'Median speedup {S["primary"]["pooled_median"]:.2f}×\n95% CI {lo:.2f}–{hi:.2f}×',transform=ax.transAxes,va='top',size=6.5)
ax.legend(handles=[Line2D([],[],color=c,marker='o',ms=3,label=l) for c,l in [(ORIG,'Original CellChat'),(ESS,'CellChatESS')]],loc='lower right',frameon=False,handlelength=1)
ax=fig.add_subplot(gs[0,1]);title(ax,'B','AutoZyme vs ESS','32 inputs · 192 paired blocks')
for st,mark,lab in zip(strata,['o','^','s'],['1k','10k-or-full','50k-or-full']):
    v=q[q.scale.eq(st)];ax.scatter(v.autozyme,v.ess,s=10,marker=mark,facecolors='none',edgecolors=ESS,lw=.5,alpha=.65,label=lab)
lim=[q[['autozyme','ess']].min().min()*.7,q[['autozyme','ess']].max().max()*1.5];ax.plot(lim,lim,'--',c=ORIG,lw=.8);ax.set_xlim(lim);ax.set_ylim(lim);log(ax,'x');log(ax)
ax.set(xlabel='AutoZyme inference (s; log)',ylabel='ESS inference time (s; log)')
ax.text(.04,.96,f'Median speedup {S["autozyme"]["pooled_median"]:.3f}×\nESS faster: {S["autozyme"]["ess_faster_blocks"]}/192',transform=ax.transAxes,va='top',size=6.3)
ax.legend(loc='lower right',frameon=False,handlelength=.8,labelspacing=.2)
ax=fig.add_subplot(gs[0,2]);title(ax,'C','Six-stage ablation','Three datasets · 3 runs/stage')
for scale,col,marker,lab in zip([1000,10000,50000],PALETTE,['o','^','s'],['1k · 10X','10k · EGAS00001004422','50k · PRJNA638151']):
    v=ab[ab.scale.eq(scale)].sort_values('stage');ax.plot(range(6),v.median_sec,marker=marker,color=col,ms=3,label=lab)
ax.set_xticks(range(6),[f'A{i}' for i in range(6)]);log(ax);ax.set(ylabel='Inference-call time (s; log)',xlabel='Ablation stage (see caption)');ax.set_ylim(.4,150000)
ax.legend(loc='upper right',bbox_to_anchor=(1,1.035),frameon=False,handlelength=1,fontsize=5.8)
groups=read('DE_groups.tsv')
network_weights=[read('Figure1D_original_weights.tsv'),read('Figure1E_ess_weights.tsv')]
max_weight=max(v.weight.max() for v in network_weights)
for i,letter,name in [(0,'D','CellChat 2.2.0.9001'),(1,'E','CellChatESS 0.2.2')]:
    ax=fig.add_subplot(gs[1,i]);title(ax,letter,name,'Glioblastoma · 1k cells · K = 9');ax.axis('off')
    pos=ax.get_position()
    netax=fig.add_axes([pos.x0-.025,pos.y0-.012,pos.width+.05,pos.height+.016])
    draw_network(netax,groups,network_weights[i],max_weight)
ax=fig.add_subplot(gs[1,2]);title(ax,'F','Numerical agreement','New 10X example · 1k cells · K = 12')
top=fr.original.max()*1.04;ax.plot([0,top],[0,top],'--',c=ORIG,lw=.7);v=fr[fr.original.ne(0)|fr.ESS.ne(0)];ax.scatter(v.original,v.ESS,s=3,c=ESS,alpha=.25,edgecolors='none');ax.set(xlim=(0,top),ylim=(0,top),xlabel='Original CellChat probability',ylabel='CellChatESS probability')
ax.set_ylabel('CellChatESS probability',labelpad=2)
ax.text(.04,.92,f'Example max |Δ| = {S["representative"]["max_abs_difference"]:.1e}',transform=ax.transAxes,size=6.2,va='top')
ax=fig.add_subplot(gs[2,0]);title(ax,'G','Peak process memory','Three methods · 192 paired blocks')
dist(ax,mem,[('Original','original',ORIG),('AutoZyme','autozyme',AUTO),('ESS','ess',ESS)],'scale',strata,sl);log(ax);ax.set(xlabel='Input stratum',ylabel='Peak RSS (GiB; log)');ax.set_ylim(.35,mem[['original','ess','autozyme']].max().max()*2)
ax.legend(handles=[Line2D([],[],color=c,marker='o',ms=3,label=l) for c,l in [(ORIG,'Original CellChat'),(AUTO,'AutoZyme'),(ESS,'ESS')]],loc='upper left',frameon=False,handlelength=.8,labelspacing=.25)
ax=fig.add_subplot(gs[2,1]);title(ax,'H','Paired RSS differences','ESS relative to AutoZyme')
dist(ax,mem,[('ESS − AutoZyme','ess_minus_autozyme_GiB',ESS)],'scale',strata,sl);ax.axhline(0,ls='--',c=AUTO,lw=.7);ax.set(xlabel='Input stratum',ylabel='ESS − AutoZyme RSS (GiB)');v=mem.ess_minus_autozyme_GiB;ax.set_ylim(v.min()*1.17,max(.04,v.max()+.04));ax.text(.02,.02,f'ESS < AutoZyme: {S["memory_ESS_lower_than_AutoZyme_blocks"]}/192 blocks',transform=ax.transAxes,size=6.3,va='bottom')
ax=fig.add_subplot(gs[2,2]);title(ax,'I','Xenium resource profile','1,156,091 cells · K = 2,357');ax.axis('off');x=S['xenium'];payload=(int(x['probability_payload_bytes'])+int(x['pvalue_payload_bytes']))/1024**3
rows=[['Inference',f'{float(x["elapsed_sec"]):,.3f} s'],['Full process',f'{x["full_process_wall_sec"]:,.2f} s'],['GNU peak RSS',f'{x["peak_RSS_GiB"]:.3f} GiB'],['Dense outputs',f'{payload:.2f} GiB']]
tab=ax.table(cellText=rows,colLabels=['Resource','Measurement'],colWidths=[.50,.50],cellLoc='left',bbox=[0,.38,1,.55]);tab.auto_set_font_size(False);tab.set_fontsize(6.9)
for (r,c),cell in tab.get_celld().items():
    cell.set_edgecolor('#cdd2d7');cell.set_linewidth(.45);cell.PAD=.04
    if r==0:cell.get_text().set_weight('bold')
    if c==1 and r>0:cell.get_text().set_ha('right')
ax.text(0,.29,f'Input expression: {int(x["input_bytes"])/1024**3:.3f} GiB\nEach dense output: {payload/2:.2f} GiB\n627 LR · 1 thread · non-spatial\nStorage per tensor: O(nLR K²)',va='top',size=7.1,linespacing=1.35)
sub=gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=gs[3,:],wspace=.10)
j=read('J_bubble_scores.tsv');order=read('J_preserved_LR_order.tsv');pairs=j[['pair','pair_rank']].drop_duplicates().sort_values('pair_rank');vmin=j.prob.min();vmax=j.prob.max()
cmap=colors.LinearSegmentedColormap.from_list('communication_blues',['#CDE2FB','#9EC5F4','#6DA7EC','#3987E5','#256ABF','#184F95','#0D366B'])
for k,engine in enumerate(['Original CellChat','CellChatESS']):
    ax=fig.add_subplot(sub[0,k]);v=j[j.engine.eq(engine)];pos=ax.get_position();ax.set_position([.15 if k==0 else .57,pos.y0,.362,pos.height]);size=(.6+(v.prob-vmin)/(vmax-vmin)*3)**2*2.2
    sc=ax.scatter(v.pair_rank,v.lr_rank,c=v.prob,s=size,cmap=cmap,vmin=vmin,vmax=vmax,lw=0);ax.set(xlim=(.25,pairs.pair_rank.max()+.75),ylim=(16.8,.2));ax.set_yticks(order.lr_rank,order.interaction_name.str.replace('_','–',regex=False) if k==0 else ['']*16);ax.set_xticks(pairs.pair_rank,pairs.pair)
    ax.tick_params(axis='both',labelsize=7.1,pad=1)
    ax.tick_params(axis='x',labelsize=6.8)
    for lab in ax.get_xticklabels():lab.set_rotation(78);lab.set_rotation_mode('anchor');lab.set_ha('right');lab.set_va('center')
    ax.grid(color='#eeeeee',lw=.35);ax.set_axisbelow(True);ax.set_title('CellChat 2.2.0.9001' if k==0 else 'CellChatESS 0.2.2',fontsize=7.5,pad=5)
    if k==0:ax.text(-.15,1.05,'J',transform=ax.transAxes,weight='bold',size=11)
    else:
        cax=ax.inset_axes([1.012,.06,.026,.80]);cb=fig.colorbar(sc,cax=cax);cb.solids.set_rasterized(False);raw=np.array([.001,.01,.05,.10,.30]);ticks=-1/np.log(raw);cb.set_ticks(ticks);cb.set_ticklabels([f'{v:g}' for v in raw]);cb.ax.tick_params(labelsize=5.8,width=.4,length=2)
        ax.text(1.125,1.015,'Communication\nprobability',transform=ax.transAxes,ha='right',va='bottom',size=6.1)
save(fig,'Figure1_final_CC22_ESS022_bold_20261101')
print('Exported updated Figure 1 into',F)

