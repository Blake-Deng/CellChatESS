#!/usr/bin/env python3
"""Portable entry point for the deposited paired supplemental figure."""
from pathlib import Path
import argparse, csv, json, shutil, sys
HERE=Path(__file__).resolve().parent
P=argparse.ArgumentParser(description=__doc__);P.add_argument('--output',type=Path,required=True);A=P.parse_args()
OUT=A.output.resolve()
if OUT.exists():P.error('Choose a new output directory')
(OUT/'figures/plot_data').mkdir(parents=True)
for f in (HERE/'response_R1_plot_data').glob('*.tsv'):shutil.copy2(f,OUT/'figures/plot_data'/f.name)
shutil.copy2(HERE/'native_code/response_R1_summary_final.json',OUT/'summary_final.json')
EVIDENCE=HERE.parent/'xenium'
with (EVIDENCE/'R2_xenium_memory_time_comparison.tsv').open() as f:resource_rows=list(csv.DictReader(f,delimiter='\t'))
source=HERE/'native_code/generated_figure_source.py'
sys.argv=[str(source),'--delivery',str(OUT)]
scope={'__name__':'__main__','__file__':str(source),'EVIDENCE':EVIDENCE,'resource_rows':resource_rows}
exec(compile(source.read_text(),str(source),'exec'),scope)
print('Rebuilt supplemental artwork in',OUT)
