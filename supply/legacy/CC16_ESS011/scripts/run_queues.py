#!/usr/bin/env python3
"""Optional independent-process reproduction. Does nothing without --execute.

Run from repository root after restoring the omitted prepared inputs and installing
version-specific R libraries. Existing evidence is never overwritten.
"""
import argparse,itertools,os,subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('mode',choices=['benchmark','compatible']);p.add_argument('--execute',action='store_true');p.add_argument('--jobs',type=int,default=8);p.add_argument('--rscript',default='Rscript');p.add_argument('--output',type=Path,default=Path('reproduction_runs'));args=p.parse_args()
if not args.execute: p.exit(message='Plan only. Add --execute after restoring inputs; no benchmark was run.\n')
out=args.output.resolve()
if out.exists():raise SystemExit('Choose a new --output directory; refusing to overwrite')
out.mkdir(parents=True);env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','RCPP_PARALLEL_NUM_THREADS']:env[k]='1'
def run(cmd,log):
 with log.open('w') as f:subprocess.run(cmd,cwd=R,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
def original_pair(row):
 id=row.experiment_id;dst=out/id;dst.mkdir()
 for method in ['original','ess']:
  run([args.rscript,str(R/'scripts/common_timer_worker.R'),'--method',method,'--id',id,'--input',str(R/'data/inputs'/id/'prepared_cellchat.rds'),'--reference',str(R/'data/references'/id/'baseline.computed.rds'),'--out',str(dst/(method+'.tsv')),'--result',str(dst/(method+'.rds'))],dst/(method+'.log'))
def compatible_input(row):
 dst=out/(row.dataset_id+'__'+row.scale_label);dst.mkdir()
 ref='-'
 for block,order in enumerate([('original','autozyme','ess'),('autozyme','ess','original'),('ess','original','autozyme'),('original','ess','autozyme'),('ess','autozyme','original'),('autozyme','original','ess')],1):
  for method in order:
   name=f'block{block}_{method}';result=dst/(name+'.rds')
   run(['/usr/bin/time','-v','-o',str(dst/(name+'.time.txt')),args.rscript,str(R/'scripts/run_compatible_method.R'),method,str(R/row.input_path),ref,str(dst/(name+'.tsv')),str(result)],dst/(name+'.log'))
   if method=='original' and ref=='-':ref=str(result)
if args.mode=='benchmark':
 data=pd.read_csv(R/'results/benchmark_216/common_external_timer_216.tsv',sep='\t')
 with ThreadPoolExecutor(max_workers=min(args.jobs,8)) as pool:list(pool.map(original_pair,data.itertuples()))
else:
 data=pd.read_csv(R/'results/autozyme_192/inputs_32.tsv',sep='\t')
 with ThreadPoolExecutor(max_workers=min(args.jobs,4)) as pool:list(pool.map(compatible_input,data.itertuples()))
