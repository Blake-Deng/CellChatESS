#!/usr/bin/env python3
"""Install the four frozen software sources into fresh, separate R libraries."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,tarfile,time
ROOT=Path(__file__).resolve().parents[1]
P=argparse.ArgumentParser();P.add_argument('--output',type=Path,required=True);P.add_argument('--dependency-library',type=Path,required=True);P.add_argument('--r',default='R');A=P.parse_args()
if A.output.exists():P.error('Choose a fresh installation output directory')
if not A.dependency_library.is_dir():P.error('Provide the locked dependency library, including CellChat/AutoZyme dependencies and expm for CellChat 1.6.1')
A.output=A.output.resolve();A.output.mkdir(parents=True);deps=A.dependency_library.resolve()
R=shutil.which(A.r)
if R is None:P.error('R not found')
Rscript=Path(R).parent/'Rscript'
if not Rscript.exists():P.error('Rscript not found beside selected R')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def extract(a,dest):
 dest.mkdir(parents=True)
 with tarfile.open(a) as tf:
  for m in tf.getmembers():
   p=(dest/m.name).resolve()
   if not p.is_relative_to(dest.resolve()):raise ValueError('Unsafe archive path')
   if m.issym() or m.islnk():
    if not (p.parent/m.linkname).resolve().is_relative_to(dest.resolve()):raise ValueError('Unsafe archive link')
  tf.extractall(dest)
 return next(p for p in dest.iterdir() if p.is_dir())
sources={}
for name in ['CellChat_2.2.0.9001','AutoZyme_0.3.1','CellChat_1.6.1']:
 j=json.loads((ROOT/'upstream'/(name+'_provenance.json')).read_text())
 archive=next((ROOT/'upstream').glob('*'+j['commit']+'.tar.gz'))
 if sha(archive)!=j['sha256']:raise SystemExit('Upstream archive hash mismatch: '+name)
 sources[name]=extract(archive,A.output/'extracted'/name)
ess=ROOT/'software/CellChatESS_0.2.2.tar.gz'
provenance=json.loads((ROOT/'audit/final_candidate.json').read_text())
if sha(ess)!=provenance['source_tarball_sha256']:raise SystemExit('ESS source tarball hash mismatch')
libs={k:A.output/'library'/k for k in ['cellchat22','autozyme','ess022','cellchat16']}
env=os.environ.copy();env['R_LIBS']=str(deps);env['AUTOZYME_DISABLED']='true'
commands=[('cellchat22',sources['CellChat_2.2.0.9001']),('autozyme',sources['AutoZyme_0.3.1']/'autozyme_r'),('ess022',ess),('cellchat16',sources['CellChat_1.6.1'])]
records=[]
for key,source in commands:
 libs[key].mkdir(parents=True)
 cmd=[R,'CMD','INSTALL','--preclean','--clean','--library='+str(libs[key]),str(source)]
 (A.output/(key+'_command.json')).write_text(json.dumps(cmd,indent=2)+'\n')
 with (A.output/(key+'_install.log')).open('w') as f:r=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
 records.append({'software':key,'exit_code':r.returncode,'library':str(libs[key])})
 (A.output/'installation_status.json').write_text(json.dumps(records,indent=2)+'\n')
 if r.returncode:raise SystemExit('Installation failed for '+key+'; inspect its log. Preserve version locks when resolving dependencies.')
rt={'rscript':str(Rscript),'ess_library':str(libs['ess022']),'cellchat_22_library':str(libs['cellchat22']),'autozyme_library':str(libs['autozyme']),'cellchat_16_library':str(libs['cellchat16']),'dependency_library':str(deps),'ess_source_id':provenance['numeric_source_id']}
(A.output/'runtime.json').write_text(json.dumps(rt,indent=2)+'\n')
print('Installed all four source packages. Next run version preflight and three-method smoke. Runtime:',A.output/'runtime.json')
