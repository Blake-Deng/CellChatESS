#!/usr/bin/env python3
"""Recompute light summaries from saved rows; no inference is run."""
from pathlib import Path
import random, statistics as st
import numpy as np,pandas as pd
from scipy.stats import spearmanr,rankdata
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/recomputed';OUT.mkdir(exist_ok=True)
def pct(v,q):return float(np.quantile(v,q,method='linear'))
a=pd.read_csv(ROOT/'results/benchmark_216/common_external_timer_216.tsv',sep='\t')
ids=sorted(a.dataset_id.unique());by={ds:a.loc[a.dataset_id==ds,'paired_speedup'].tolist() for ds in ids};med={ds:st.median(by[ds]) for ds in ids}
rng=random.Random(20261026);pooled=[];dataset=[]
for _ in range(10000):
 sampled=[rng.choice(ids) for _ in ids];pooled.append(st.median(x for ds in sampled for x in by[ds]));dataset.append(st.median(med[ds] for ds in sampled))
z=pd.DataFrame([dict(estimand='median_of_216_paired_inference_call_ratios',estimate=st.median(a.paired_speedup),ci95_lower=pct(pooled,.025),ci95_upper=pct(pooled,.975)),dict(estimand='median_of_12_dataset_medians',estimate=st.median(med.values()),ci95_lower=pct(dataset,.025),ci95_upper=pct(dataset,.975))]);z.to_csv(OUT/'bootstrap_216.tsv',sep='\t',index=False)
ref=pd.read_csv(ROOT/'results/benchmark_216/common_external_timer_bootstrap_summary.tsv',sep='\t');assert np.allclose(z[['estimate','ci95_lower','ci95_upper']],ref[['estimate','ci95_lower','ci95_upper']])
d=pd.read_csv(ROOT/'results/dataset_characteristics/dataset_characteristics_216.tsv',sep='\t');rows=[]
for col in ['matrix_nonzero','matrix_elements','data_signaling_size_mib','lr_group_pairs_after_pruning','f_LR']:rows.append(dict(analysis='Spearman 216 paired records',feature=col,n=216,rho=spearmanr(d[col],d.paired_speedup).statistic))
g=d.groupby('dataset_id')[['f_LR','paired_speedup']].median();rows.append(dict(analysis='Spearman 12 dataset medians',feature='f_LR',n=12,rho=spearmanr(g.f_LR,g.paired_speedup).statistic))
x,y,z=rankdata(d.f_LR),rankdata(d.paired_speedup),rankdata(d.matrix_nonzero);X=np.column_stack([np.ones(len(d)),z]);rx=x-X@np.linalg.lstsq(X,x,rcond=None)[0];ry=y-X@np.linalg.lstsq(X,y,rcond=None)[0];rows.append(dict(analysis='Partial rank; control matrix_nonzero',feature='f_LR',n=216,rho=np.corrcoef(rx,ry)[0,1]));pd.DataFrame(rows).to_csv(OUT/'correlations.tsv',sep='\t',index=False)
print('Bootstrap and correlations reproduced from existing records.')
