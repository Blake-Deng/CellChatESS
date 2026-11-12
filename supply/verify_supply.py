#!/usr/bin/env python3
"""Verify deposited files and unchanged algorithm source; no R execution needed."""
from pathlib import Path
import argparse, csv, hashlib, json, subprocess, sys

SUP=Path(__file__).resolve().parent
ROOT=SUP.parent
P=argparse.ArgumentParser(description=__doc__)
P.add_argument('--refresh-manifest',action='store_true',help='Explicitly regenerate after an intentional authoring change')
A=P.parse_args()

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))

manifest=SUP/'provenance/supply_SHA256_manifest.tsv'
ignored={manifest,SUP/'provenance/verification_report.json'}
files=sorted(p for p in SUP.rglob('*') if p.is_file() and p not in ignored
             and '__pycache__' not in p.parts and p.name!='.DS_Store')
if A.refresh_manifest:
    with manifest.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['path','size_bytes','sha256'],delimiter='\t')
        w.writeheader()
        for p in files:w.writerow({'path':str(p.relative_to(ROOT)), 'size_bytes':p.stat().st_size,'sha256':sha(p)})
if not manifest.exists():raise SystemExit('Manifest missing; author must generate it after assembling the final document.')
expected=rows(manifest)
errors=[]
for r in expected:
    p=ROOT/r['path']
    if not p.is_file():errors.append('Missing: '+r['path'])
    elif p.stat().st_size!=int(r['size_bytes']) or sha(p)!=r['sha256']:errors.append('Hash/size mismatch: '+r['path'])
if {r['path'] for r in expected}!={str(p.relative_to(ROOT)) for p in files}:errors.append('Supply file inventory differs from manifest')

source=rows(SUP/'provenance/candidate_source_manifest.tsv')
documentary={'README.md','.Rbuildignore','CITATION.cff','.zenodo.json','NEWS.md'}
for r in source:
    p=SUP/'provenance/baseline_documentation'/r['path'] if r['path'] in documentary else ROOT/r['path']
    if not p.is_file() or sha(p)!=r['sha256']:errors.append('Baseline source mismatch: '+r['path'])
if '^supply$' not in (ROOT/'.Rbuildignore').read_text().splitlines():errors.append('Missing R build exclusion for supply/')
if 'Version: 0.2.2' not in (ROOT/'DESCRIPTION').read_text():errors.append('Wrong package version')

refs=rows(SUP/'provenance/result_version_build_map.tsv')
for r in refs:
    if not (ROOT/r['evidence']).exists():errors.append('Missing mapped evidence: '+r['evidence'])
required={'input_characteristics_216.tsv':216,'prepared_input_seeds.tsv':216,
          'primary_pairs_216.tsv':216,'autozyme_pairs_192.tsv':192,'autozyme_inputs_32.tsv':32,
          'all_main_measurements_1008.tsv':1008}
for path,n in required.items():
    actual=len(rows(SUP/'data'/path))
    if actual!=n:errors.append(f'Row count: {path}: {actual} != {n}')
# Historical AutoZyme evidence stays separate and keeps its recorded versions.
legacy=SUP/'legacy/CC16_AutoZyme_original_submission'
for name in ['autozyme_216_all.tsv','ess_216_rerun_all.tsv','autozyme_vs_ess_peak_216.tsv']:
    path=legacy/name
    if not path.is_file():
        errors.append('Missing historical comparison: '+name)
        continue
    historical=rows(path)
    if len(historical)!=216:errors.append('Historical comparison row count: '+name)
    for row in historical:
        if name=='autozyme_vs_ess_peak_216.tsv':
            valid=row['package_version_az']=='0.3.1' and row['package_version_ess']=='0.1.1'
        else:
            valid=row['CellChat_version']=='1.6.1' and row['CellChatESS_version']=='0.1.1' and row['autozyme_version']=='0.3.1'
        if not valid:
            errors.append('Historical comparison version mismatch: '+name)
            break
for forbidden in ['S_H_Xenium_resources.tsv','S_H2_stage_events.tsv','S_H3_RSS_samples.tsv','S_I_result_version_build_map.tsv']:
    if (SUP/'data'/forbidden).exists():errors.append('Superseded table in current data/: '+forbidden)
for p in files:
    if p.suffix.lower() in {'.so','.o','.dll','.dylib','.rds'}:errors.append('Excluded binary/raw R object: '+str(p.relative_to(ROOT)))
report={'status':'PASS' if not errors else 'FAIL','supply_files':len(expected),
        'supply_bytes':sum(int(r['size_bytes']) for r in expected),'baseline_source_files_verified':len(source),
        'unchanged_computational_source':True if not any('Baseline source' in e for e in errors) else False,
        'baseline_commit':'139a61c33d4c09568f48f4afd0a16dc57bb22f4f',
        'baseline_tree':'62fb33d233e8b2296eca11b7eefade7745429e10',
        'documentary_changes_allowed':sorted(documentary),'errors':errors,
        'note':'Verifies deposited records; does not rerun the scientific computations or recover excluded large RDS outputs.'}
if (ROOT/'.git').is_dir():
    try:subprocess.run(['git','-C',str(ROOT),'cat-file','-e',report['baseline_commit']+'^{commit}'],check=True,capture_output=True)
    except subprocess.CalledProcessError:errors.append('Baseline commit missing from local Git history')
(SUP/'provenance/verification_report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
sys.exit(bool(errors))
