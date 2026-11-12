"""Update only R1D in the existing native scientific plotting source."""
from pathlib import Path
import csv
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parent.parent
WORK = Path(__file__).resolve().parent
OLD = ROOT / 'Response_FigureR1_bold_20261101'
EVIDENCE = ROOT / 'CellChatESS_0.2.2_CellChat_2.2.0.9001/Xenium_627LR_evidence_minimal_20261102/Xenium_627LR_evidence'
DELIVERY = WORK / 'figure'
DATA = DELIVERY / 'figures/plot_data'
DATA.mkdir(parents=True, exist_ok=True)

for name in ['A_primary_pairs_216.tsv', 'B_autozyme_pairs_192.tsv', 'R1C_high_retention_pairs.tsv']:
    (DATA / name).write_bytes((OLD / 'figures/plot_data' / name).read_bytes())
summary = json.loads((OLD / 'summary_final.json').read_text())
resource_rows = list(csv.DictReader((EVIDENCE / 'R2_xenium_memory_time_comparison.tsv').open(), delimiter='\t'))
summary['xenium'] = {'scope': 'Full 627-LR paired runs; clean performance builds', 'methods': resource_rows}
(DELIVERY / 'summary_final.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n')

source = (OLD / 'redraw_response_R1.py').read_text()
source = source[:source.index('sub=gs[1,1].subgridspec')]
source = source.replace("'data_changes':False", "'data_changes':True,'changed_panels':['D']")
source = source.replace('Pinned CellChat 2.2.0.9001 / ESS 0.2.2 results; typography and layout update only',
                        'Panels A-C unchanged; panel D uses complete paired Xenium memory records')
source += '''
sub=gs[1,1].subgridspec(2,1,height_ratios=[1,1],hspace=.48)
paired_sources=[]
for i,(method,label,color) in enumerate([
    ('original','Original CellChat',ORIG),('ess','ESS',ESS)]):
    folder=EVIDENCE/method
    raw=pd.read_csv(folder/'rss_samples.tsv',sep='\\t',float_precision='round_trip')
    stages=pd.read_csv(folder/'stage_events.tsv',sep='\\t',float_precision='round_trip')
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
    raw.to_csv(D/f'R1D_{method}_paired_RSS.tsv',sep='\\t',index=False)
    stages.to_csv(D/f'R1D_{method}_paired_stages.tsv',sep='\\t',index=False)
    paired_sources.append({'method':method,'input_sha256':run['input']['input_sha256'],
        'inference_sec':run['inference_function_elapsed_sec'],
        'aggregate_sec':run['aggregateNet_elapsed_sec'],
        'GNU_peak_RSS_GiB':peak,
        'inference_end_min':infer_end,'aggregate_end_min':agg_end})
save(fig,'FigureR1_paired_CC22_20261102')
(Q/'paired_sources.json').write_text(json.dumps(paired_sources,indent=2)+'\\n')
'''
(WORK / 'generated_figure_source.py').write_text(source)
sys.argv = ['generated_figure_source.py', '--delivery', str(DELIVERY)]
exec(compile(source, str(WORK / 'generated_figure_source.py'), 'exec'),
     {'__name__': '__main__', '__file__': str(WORK / 'generated_figure_source.py'),
      'EVIDENCE': EVIDENCE, 'resource_rows': resource_rows})

# Retain numerical input identities for the three unchanged panels.
audit = {name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
         for name in ['A_primary_pairs_216.tsv', 'B_autozyme_pairs_192.tsv', 'R1C_high_retention_pairs.tsv']}
(WORK / 'unchanged_panel_data_sha256.json').write_text(json.dumps(audit, indent=2) + '\n')
