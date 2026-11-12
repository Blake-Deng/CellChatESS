"""Shared, explicit IO helpers for server follow-ups; no inference on import."""
from pathlib import Path
import csv, hashlib, json, os, re, subprocess, time
ROOT=Path(__file__).resolve().parents[1]

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def rows(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))

def table(p,data):
    data=list(data);p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    fields=list(dict.fromkeys(k for r in data for k in r))
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(data)

def runtime():return json.loads((ROOT/'server_install/runtime.json').read_text())
def env1(rt):
    e=os.environ.copy();e['R_LIBS']=rt['dependency_library'];e['AUTOZYME_DISABLED']='true'
    for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','RCPP_PARALLEL_NUM_THREADS','AUTOZYME_THREADS','ZYME_THREADS']:e[k]='1'
    e['LC_ALL']='C';return e

def resources(p):
    s=Path(p).read_text();rss=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',s)
    wall=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)',s)
    exitcode=re.search(r'Exit status:\s*(\d+)',s)
    sec=0.
    if wall:
        for x in wall[1].split(':'):sec=sec*60+float(x)
    return {'peak_RSS_GiB':int(rss[1])/1024**2 if rss else None,
            'full_process_wall_sec':sec if wall else None,'exit_status':int(exitcode[1]) if exitcode else None}

def record_run(cmd,directory,env,gnu=True):
    d=Path(directory)
    if d.exists():raise RuntimeError('Refusing to overwrite '+str(d))
    d.mkdir(parents=True)
    cmd=list(map(str,cmd));full=['/usr/bin/time','-v','-o',str(d/'gnu_time.txt')]+cmd if gnu else cmd
    (d/'command.json').write_text(json.dumps(full,indent=2)+'\n')
    start=time.time()
    with (d/'run.log').open('w') as f:rc=subprocess.run(full,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
    (d/'exit.json').write_text(json.dumps({'exit_code':rc,'started_epoch':start,'finished_epoch':time.time()},indent=2)+'\n')
    return rc

def compare(ref,other,directory,rt):
    d=Path(directory);d.mkdir(parents=True,exist_ok=True)
    cmd=[rt['rscript'],'--vanilla',str(ROOT/'reproduction/compare_results.R'),str(ref),str(other),str(d/'comparison.tsv')]
    with (d/'compare.log').open('w') as f:r=subprocess.run(cmd,cwd=ROOT,env=env1(rt),stdout=f,stderr=subprocess.STDOUT)
    return r.returncode==0
