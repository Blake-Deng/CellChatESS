#!/usr/bin/env python3
"""Finish delivery indexing after the existing strict queue exits; no inference."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]

def read_json(path):
    return json.loads(path.read_text()) if path.is_file() else {}

def write_json(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
    temporary.replace(path)

def pipeline_alive(pid):
    try:
        command = (Path('/proc') / str(pid) / 'cmdline').read_bytes()
        return b'execute_xenium_strict_pipeline.py' in command
    except OSError:
        return False

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    out = args.output.resolve()
    state_path = out / 'background_queue_status.json'
    lock = out / '.background_queue.lock'
    descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    os.write(descriptor, str(os.getpid()).encode())
    os.close(descriptor)
    started = time.time()
    pid = read_json(out / 'pipeline_process.json')['pid']
    def record(status, **values):
        state = {'status': status, 'queue_finalizer_pid': os.getpid(),
                 'existing_pipeline_pid': pid, 'started_epoch': started,
                 'updated_epoch': time.time(),
                 'does_not_launch_duplicate_inference': True, **values}
        write_json(state_path, state)
        print(json.dumps(state, ensure_ascii=False), flush=True)
    record('WAITING_FOR_EXISTING_PIPELINE')
    while True:
        state = read_json(out / 'pipeline_status.json')
        alive = pipeline_alive(pid)
        if state.get('status') != 'RUNNING' and not alive:
            break
        if not alive and state.get('status') == 'RUNNING':
            reason = 'Pipeline exited without a terminal status; no automatic rerun or completed-result claim.'
            record('FAILED_ABNORMAL_PIPELINE_EXIT', reason=reason)
            (out / 'QUEUE_FAILURE_CN.md').write_text(
                '# 后台队列异常退出\n\n主队列进程已退出，但未写入终止状态。'
                '无法确认全量完成；不能报告完整配对性能比或全量一致性。'
                '请保留并核查 pipeline.log、各方法 run.log、GNU time、RSS 和阶段记录。\n')
            return 2
        time.sleep(45)
    record('INDEXING_TERMINAL_EVIDENCE', pipeline_terminal_status=state)
    command = [sys.executable, str(ROOT / 'reproduction/audit_xenium_paired_delivery.py'),
               '--output', str(out)]
    rc = subprocess.run(command, cwd=ROOT).returncode
    integrity = read_json(out / 'terminal_delivery_integrity.json')
    expected = ['Xenium_627LR_paired_report_CN.md',
                'R2_xenium_memory_time_comparison.tsv', 'validation_acceptance.tsv',
                'xenium_paired_RSS.pdf', 'xenium_paired_RSS.svg', 'xenium_paired_RSS.png']
    files_present = {name: (out / name).is_file() for name in expected}
    complete = rc == 0 and integrity.get('full_comparison_completed', False) and all(files_present.values())
    record('COMPLETED' if complete else 'INCOMPLETE',
           audit_exit_code=rc, pipeline_terminal_status=state,
           full_pair_and_comparison_completed=integrity.get('full_comparison_completed', False),
           expected_delivery_files=files_present,
           interpretation='Completion and numerical agreement are separate; inspect validation_acceptance.tsv and comparison outputs.')
    return 0 if complete else 2

if __name__ == '__main__':
    sys.exit(main())
