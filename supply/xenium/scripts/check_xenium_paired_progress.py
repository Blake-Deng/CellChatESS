#!/usr/bin/env python3
from pathlib import Path
import argparse,json,os,time
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output
state=out/'pipeline_status.json'
print('Pipeline:',json.loads(state.read_text()) if state.exists() else 'status missing')
for method in ('original','ess'):
 d=out/method;status=d/'run_status.json';pidfile=d/'actual_process.json'
 print(method+':',json.loads(status.read_text()) if status.exists() else 'not completed')
 if pidfile.exists():
  pid=json.loads(pidfile.read_text())['R_pid'];proc=Path('/proc')/str(pid)
  if proc.exists():
   try:
    args=proc.joinpath('cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
    if 'run_xenium_paired_worker.R' in args:
     stat=proc.joinpath('stat').read_text().split();hz=os.sysconf('SC_CLK_TCK')
     elapsed=float(Path('/proc/uptime').read_text().split()[0])-int(stat[21])/hz
     cpu=(int(stat[13])+int(stat[14]))/hz
     print('  active R pid=%d elapsed=%.1f min average CPU=%.1f%%'%(pid,elapsed/60,100*cpu/elapsed))
     print('\n'.join('  '+x for x in proc.joinpath('status').read_text().splitlines() if x.startswith(('State:','VmRSS:','VmHWM:','Threads:'))))
   except (OSError,ValueError,IndexError):print('  process exited during observation')
mem={x.split(':')[0]:int(x.split()[1]) for x in Path('/proc/meminfo').read_text().splitlines()}
print('MemAvailable GiB: %.2f'%(mem['MemAvailable']/1024**2))
log=out/'pipeline.log'
if log.exists():print('\n'.join(log.read_text().splitlines()[-4:]))
for rel in ['comparison_status.json','R2_xenium_memory_time_comparison.tsv']:
 f=out/rel
 if f.exists():print(rel+'\n'+f.read_text())
