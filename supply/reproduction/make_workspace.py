#!/usr/bin/env python3
"""Restore the frozen server script layout in a fresh external workspace."""
from pathlib import Path
import argparse, csv, hashlib, json, shutil, tarfile, urllib.request

SUP=Path(__file__).resolve().parents[1]
REPO=SUP.parent
P=argparse.ArgumentParser(description=__doc__)
P.add_argument('--output',type=Path,required=True)
P.add_argument('--fetch-upstream',action='store_true',help='Download exact committed sources, with mandatory SHA256 checks')
A=P.parse_args()
OUT=A.output.resolve()
if OUT.exists():P.error('Choose a new empty output path')
if OUT==REPO or REPO in OUT.parents:P.error('Place the workspace outside the source repository')
OUT.mkdir(parents=True)
for name in ['reproduction','validation']:
    shutil.copytree(SUP/'reproduction/frozen'/name,OUT/name)
src=OUT/'software/CellChatESS';src.mkdir(parents=True)
with (SUP/'provenance/candidate_source_manifest.tsv').open() as f:
    records=list(csv.DictReader(f,delimiter='\t'))
for r in records:
    p=(SUP/'provenance/baseline_documentation'/r['path']
       if r['path'] in {'README.md','.Rbuildignore','CITATION.cff','.zenodo.json'} else REPO/r['path'])
    if hashlib.sha256(p.read_bytes()).hexdigest()!=r['sha256']:raise SystemExit('Source hash mismatch: '+r['path'])
    dest=src/r['path'];dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
shutil.copy2(SUP/'provenance/archives/CellChatESS_0.2.2_source.tar.gz',OUT/'software/CellChatESS_0.2.2.tar.gz')
for name in ['audit','upstream','environment','logs','results/new']:(OUT/name).mkdir(parents=True,exist_ok=True)
for name in ['candidate_source_manifest.tsv','final_candidate.json','server_source_identity.json']:
    shutil.copy2(SUP/'provenance'/name,OUT/'audit'/name)
shutil.copytree(SUP/'environment/server',OUT/'environment/server')
for p in (SUP/'provenance/upstream').glob('*.json'):
    shutil.copy2(p,OUT/'upstream'/p.name)
    j=json.loads(p.read_text())
    if A.fetch_upstream:
        url=j['archive_url'];dest=OUT/'upstream'/url.rsplit('/',1)[-1]
        # Frozen installers locate by commit + .tar.gz, independent of URL basename.
        dest=OUT/'upstream'/(j['repo'].split('/')[-1]+'-'+j['commit']+'.tar.gz')
        urllib.request.urlretrieve(url,dest)
        if hashlib.sha256(dest.read_bytes()).hexdigest()!=j['sha256']:raise SystemExit('Upstream SHA256 mismatch: '+p.name)
print('Created frozen-script workspace:',OUT)
print('Edit reproduction/config/runtime.example.json and use an explicit input path map.')
print('Large prepared inputs and package dependencies are not downloaded by this command.')
