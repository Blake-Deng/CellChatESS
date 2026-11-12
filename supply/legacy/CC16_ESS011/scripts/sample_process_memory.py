#!/usr/bin/env python3
"""Attach /proc RSS/PSS sampling to the profiled R process already in flight."""
import csv
import datetime as dt
import time
from pathlib import Path
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('--pid',type=int,required=True)
parser.add_argument('--stages',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
STAGES=args.stages;OUT=args.output
if OUT.exists():raise SystemExit('Refusing to overwrite existing RSS samples')
OUT.parent.mkdir(parents=True,exist_ok=True)
SCRIPT_TOKEN = 'run_xenium_profiled.R'

def find_r_pid():
    for entry in Path('/proc').iterdir():
        if not entry.name.isdigit():
            continue
        try:
            cmdline = (entry / 'cmdline').read_bytes().replace(b'\0', b' ').decode(errors='replace')
            comm = (entry / 'comm').read_text().strip()
        except OSError:
            continue
        if SCRIPT_TOKEN in cmdline and comm in ('R', 'Rscript'):
            return int(entry.name)
    return None

def latest_stage():
    latest = 'process_startup'
    try:
        with STAGES.open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                if row.get('stage'):
                    latest = row['stage']
    except OSError:
        pass
    return latest

def read_status(pid, names):
    vals = {name: '' for name in names}
    try:
        for line in Path(f'/proc/{pid}/status').read_text().splitlines():
            for name in names:
                if line.startswith(name + ':'):
                    vals[name] = line.split(':', 1)[1].strip().split()[0]
    except OSError:
        pass
    return vals

def read_smaps(pid):
    names = ['Rss', 'Pss', 'Pss_Anon', 'Pss_File', 'Private_Clean', 'Private_Dirty',
             'Shared_Clean', 'Shared_Dirty', 'Anonymous', 'Swap']
    vals = {name: '' for name in names}
    try:
        for line in Path(f'/proc/{pid}/smaps_rollup').read_text().splitlines():
            for name in names:
                if line.startswith(name + ':'):
                    vals[name] = line.split(':', 1)[1].strip().split()[0]
    except OSError:
        pass
    return vals

fields = ['timestamp_utc', 'elapsed_sec', 'pid', 'stage', 'VmRSS_kB', 'RssAnon_kB', 'RssFile_kB',
          'VmSwap_kB', 'smaps_Rss_kB', 'smaps_Pss_kB', 'Pss_Anon_kB', 'Pss_File_kB',
          'Private_Clean_kB', 'Private_Dirty_kB', 'Shared_Clean_kB', 'Shared_Dirty_kB',
          'Anonymous_kB', 'smaps_Swap_kB']
start = time.monotonic()
with OUT.open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    handle.flush()
    pid = args.pid
    while pid is None:
        pid = find_r_pid()
        if pid is None:
            time.sleep(0.5)
    while Path(f'/proc/{pid}').exists():
        status = read_status(pid, ['VmRSS', 'RssAnon', 'RssFile', 'VmSwap'])
        smaps = read_smaps(pid)
        record = {name: '' for name in fields}
        record.update({
            'timestamp_utc': dt.datetime.now(dt.timezone.utc).isoformat(timespec='milliseconds'),
            'elapsed_sec': f'{time.monotonic()-start:.3f}', 'pid': pid, 'stage': latest_stage(),
            'VmRSS_kB': status['VmRSS'], 'RssAnon_kB': status['RssAnon'],
            'RssFile_kB': status['RssFile'], 'VmSwap_kB': status['VmSwap'],
            'smaps_Rss_kB': smaps['Rss'], 'smaps_Pss_kB': smaps['Pss'],
            'Pss_Anon_kB': smaps['Pss_Anon'], 'Pss_File_kB': smaps['Pss_File'],
            'Private_Clean_kB': smaps['Private_Clean'], 'Private_Dirty_kB': smaps['Private_Dirty'],
            'Shared_Clean_kB': smaps['Shared_Clean'], 'Shared_Dirty_kB': smaps['Shared_Dirty'],
            'Anonymous_kB': smaps['Anonymous'], 'smaps_Swap_kB': smaps['Swap'],
        })
        writer.writerow(record)
        handle.flush()
        time.sleep(1.0)
print(f'RSS sampling ended for PID {pid}; output={OUT}')
