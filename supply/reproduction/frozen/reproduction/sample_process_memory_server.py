#!/usr/bin/env python3
"""Bind a VmRSS trace to a single R PID and record actual acquisition times."""
from pathlib import Path
import argparse,csv,time,datetime,json
p=argparse.ArgumentParser();p.add_argument('--pid',type=int,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--stages',type=Path,required=True);p.add_argument('--interval',type=float,default=1.0);a=p.parse_args()
if a.output.exists():raise SystemExit('Fresh sample output required')
proc=Path('/proc')/str(a.pid);start_id=proc.joinpath('stat').read_text().split()[21]
fields=['timestamp_utc','sample_start_epoch_ms','sample_end_epoch_ms','elapsed_from_sampler_start_s','pid','stage','VmRSS_kB','VmHWM_kB','VmSwap_kB','sampling_read_duration_s']
start=time.monotonic();last_stage='process_startup'
with a.output.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();f.flush()
    while proc.exists():
        t0=time.time();m0=time.monotonic()
        try:
            st=proc.joinpath('stat').read_text().split()
            if st[21]!=start_id or st[2]=='Z':break
            status={line.split(':',1)[0]:line.split(':',1)[1].strip().split()[0] for line in proc.joinpath('status').read_text().splitlines() if ':' in line and line.split(':',1)[1].strip()}
        except (OSError,IndexError):break
        try:
            for r in csv.DictReader(a.stages.open(),delimiter='\t'):
                if r.get('stage'):last_stage=r['stage']
        except OSError:pass
        t1=time.time()
        w.writerow(dict(timestamp_utc=datetime.datetime.fromtimestamp(t0,datetime.timezone.utc).isoformat(timespec='milliseconds'),sample_start_epoch_ms=f'{t0*1000:.3f}',sample_end_epoch_ms=f'{t1*1000:.3f}',elapsed_from_sampler_start_s=f'{time.monotonic()-start:.6f}',pid=a.pid,stage=last_stage,VmRSS_kB=status.get('VmRSS',''),VmHWM_kB=status.get('VmHWM',''),VmSwap_kB=status.get('VmSwap',''),sampling_read_duration_s=f'{time.monotonic()-m0:.6f}'));f.flush()
        time.sleep(max(0,a.interval-(time.monotonic()-m0)))
a.output.with_suffix('.identity.json').write_text(json.dumps({'pid':a.pid,'proc_start_ticks':start_id,'requested_interval_sec':a.interval,'source':'/proc/<PID>/status VmRSS; no smoothing or interpolation','clock':'Unix epoch sampling start and end; align separately to stage event'},indent=2)+'\n')
