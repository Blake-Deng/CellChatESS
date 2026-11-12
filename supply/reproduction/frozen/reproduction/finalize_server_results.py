#!/usr/bin/env python3
"""Summarize completed CC22 timing queues and extract new feature records."""
from pathlib import Path
import csv,json,subprocess,os
ROOT=Path(__file__).resolve().parents[1]
rt=json.loads((ROOT/'server_install/runtime.json').read_text())
out=ROOT/'results/new/core22_statistics'
if out.exists(): raise SystemExit('Refusing to overwrite statistics output')
merged=ROOT/'results/new/primary216_merged'
primary=merged if merged.exists() else ROOT/'results/new/primary216'
auto=ROOT/'results/new/autozyme192'
subprocess.run(['python3',str(ROOT/'reproduction/summarize_and_plot.py'),'--primary',str(primary),'--autozyme',str(auto),'--output',str(out)],check=True)
feat=out/'dataset_characteristics_216.tsv'
env={**os.environ,'R_LIBS':rt['dependency_library'],'AUTOZYME_DISABLED':'true'}
subprocess.run([rt['rscript'],'--vanilla',str(ROOT/'reproduction/extract_dataset_characteristics.R'),str(ROOT/'reproduction/config/primary_216_CC22.tsv'),str(feat),rt['ess_library'],rt['cellchat_22_library']],env=env,check=True)
import pandas as pd
f=pd.read_csv(feat,sep='\t');t=pd.read_csv(out/'plot_data/primary_pairs.tsv',sep='\t')
f=f.merge(t[['experiment_id','original_over_ess']],on='experiment_id',how='left',validate='one_to_one');f.to_csv(feat,sep='\t',index=False)
(out/'feature_definitions.txt').write_text('matrix_nonzero is the count of nonzero entries in data.signaling. f_LR = sum_l |{(s,r): L_l(s)>0 and R_adj_l(r)>0}|/(nLR*K^2). f_union is the count of distinct supported group pairs divided by K^2.\n')
print('NEW_CORE_STATISTICS_COMPLETE',out)
