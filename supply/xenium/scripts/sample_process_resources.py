#!/usr/bin/env python3
"""One-second /proc sampler bound to a single R process start-time identity."""
from __future__ import annotations
import argparse, csv, datetime, json, time
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument('--pid',type=int,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--stages',type=Path,required=True)
p.add_argument('--interval',type=float,default=1.0)
a=p.parse_args()
proc=Path('/proc')/str(a.pid)
if a.output.exists(): raise SystemExit(f'refusing to overwrite {a.output}')
start_ticks=proc.joinpath('stat').read_text().split()[21]
fields=['timestamp_utc','sample_start_epoch_ms','sample_end_epoch_ms','elapsed_from_sampler_start_s',
        'pid','stage','VmRSS_kB','VmHWM_kB','VmSwap_kB','Threads','sampling_read_duration_s']
start=time.monotonic(); last_stage='process_startup'
with a.output.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();f.flush()
    while proc.exists():
        t0=time.time();m0=time.monotonic()
        try:
            st=proc.joinpath('stat').read_text().split()
            if st[21]!=start_ticks or st[2]=='Z': break
            status={line.split(':',1)[0]:line.split(':',1)[1].strip().split()[0]
                    for line in proc.joinpath('status').read_text().splitlines()
                    if ':' in line and line.split(':',1)[1].strip()}
        except (OSError,IndexError): break
        try:
            with a.stages.open() as sf:
                for row in csv.DictReader(sf,delimiter='\t'):
                    if row.get('stage'):last_stage=row['stage']
        except OSError:pass
        t1=time.time()
        w.writerow({'timestamp_utc':datetime.datetime.fromtimestamp(t0,datetime.timezone.utc).isoformat(timespec='milliseconds'),
          'sample_start_epoch_ms':f'{t0*1000:.3f}','sample_end_epoch_ms':f'{t1*1000:.3f}',
          'elapsed_from_sampler_start_s':f'{time.monotonic()-start:.6f}','pid':a.pid,'stage':last_stage,
          'VmRSS_kB':status.get('VmRSS',''),'VmHWM_kB':status.get('VmHWM',''),
          'VmSwap_kB':status.get('VmSwap',''),'Threads':status.get('Threads',''),
          'sampling_read_duration_s':f'{time.monotonic()-m0:.6f}'})
        f.flush()
        time.sleep(max(0,a.interval-(time.monotonic()-m0)))
a.output.with_suffix('.identity.json').write_text(json.dumps({'pid':a.pid,'proc_start_ticks':start_ticks,
 'requested_interval_sec':a.interval,'source':'/proc/<PID>/status; no smoothing or interpolation',
 'clocks':'Unix epoch start/end timestamps plus monotonic elapsed time'},indent=2)+'\n')
