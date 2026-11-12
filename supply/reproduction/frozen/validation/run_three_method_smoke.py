#!/usr/bin/env python3
"""Run each implementation in a fresh R process; no benchmark claims."""
from pathlib import Path
import argparse,subprocess,json,hashlib,os
P=argparse.ArgumentParser();P.add_argument('--runtime',type=Path,required=True);P.add_argument('--input',type=Path,required=True);P.add_argument('--output',type=Path,required=True);P.add_argument('--raw-use',choices=['TRUE','FALSE'],default='TRUE');A=P.parse_args()
ROOT=Path(__file__).resolve().parents[1];rt=json.loads(A.runtime.read_text())
if A.output.exists():P.error('Choose a fresh output directory')
A.output.mkdir(parents=True);A.input=A.input.resolve()
(A.output/'input_sha256.txt').write_text(hashlib.sha256(A.input.read_bytes()).hexdigest()+'  '+str(A.input)+'\n')
for method in ['original','autozyme','ess','dense']:
 out=A.output/method;out.mkdir()
 cmd=[rt['rscript'],'--vanilla',str(ROOT/'validation/smoke_worker.R'),method,str(A.input),str(out.resolve()),rt['ess_library'],rt['cellchat_22_library'],rt['autozyme_library'],rt.get('dependency_library',''),'2.2.0.9001',A.raw_use]
 (out/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
 with (out/'run.log').open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 if r.returncode:raise SystemExit('FAILED '+method+'; inspect '+str(out/'run.log'))
 print('DONE',method,flush=True)
r=subprocess.run([rt['rscript'],'--vanilla',str(ROOT/'validation/compare_smoke.R'),str(A.output.resolve())])
raise SystemExit(r.returncode)
