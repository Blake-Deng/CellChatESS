#!/usr/bin/env python3
"""Index original time logs without changing them or running any experiments."""
from pathlib import Path
import csv
import re
import hashlib

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    with (ROOT / path).open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def main():
    rows = []
    for row in read('results/autozyme_192/run_manifest_576.tsv'):
        raw = read(row['raw_tsv'])[0]
        eid = Path(raw['input_path']).parent.name
        rows.append({
            'log_path': row['time_path'], 'experiment_id': row['run_id'],
            'dataset_id': row['dataset_id'], 'input_id': eid,
            'input_sha256': row['input_sha256'], 'method': row['method'],
            'CellChat': row['CellChat'],
            'CellChatESS': row['ESS'] if row['method'] == 'ess' else 'not used',
            'AutoZyme': row['AutoZyme'] if row['method'] == 'autozyme' else 'not used',
            'environment': 'compatible: CellChat 2.2.0.9001 + ESS 0.2.2 + AutoZyme 0.3.1',
            'record_path': row['raw_tsv'], 'note': 'package versions in method context; other installed methods not used',
        })
    ablation = read('results/ablation/six_stage_ablation_18.tsv')
    for path in sorted((ROOT / 'results/ablation/raw_process_logs').glob('*.time.txt')):
        scale, method = re.fullmatch(r'rss_(1k|10k|50k)_(.+)\.time\.txt', path.name).groups()
        matches = [r for r in ablation if r['scale'] == scale and r['method'] == method]
        assert len(matches) == 1
        row = matches[0]
        stage = row['stage_code']
        version = 'custom 0.1.1.9000 (archived intermediate; not released)' if stage in ['A3', 'A4'] else '0.1.1'
        note = f'{stage}: ' + ('archived intermediate driver; not a released implementation' if stage == 'A2' else row['implementation'])
        rows.append({
            'log_path': str(path.relative_to(ROOT)), 'experiment_id': row['experiment_id'],
            'dataset_id': row['experiment_id'].split('__cells-')[0],
            'input_id': row['experiment_id'], 'input_sha256': 'not recorded in this process log',
            'method': method, 'CellChat': '1.6.1', 'CellChatESS': version,
            'AutoZyme': 'not used', 'environment': 'main; archived ablation variants explicitly labelled',
            'record_path': 'results/ablation/raw_records/' + path.name.replace('.time.txt', '.tsv'),
            'note': note,
        })
    rows.append({
        'log_path': 'results/xenium/GNU_time.txt',
        'experiment_id': 'Xenium_grid50_ESS022_resource_profile',
        'dataset_id': 'Xenium', 'input_id': 'Xenium_grid50_1156091cells_K2357_nLR627',
        'input_sha256': 'see results/LARGE_FILES.tsv; no input hash in GNU time record',
        'method': 'ess_profiled', 'CellChat': '2.2.0.9001',
        'CellChatESS': '0.2.2 (instrumented profiling build)', 'AutoZyme': 'not used',
        'environment': 'compatible; instrumented released ESS 0.2.2 for observation',
        'record_path': 'results/xenium/run_summary.tsv',
        'note': 'full-process wall time/RSS, not inference-call timer',
    })
    logs = {str(p.relative_to(ROOT)) for p in ROOT.rglob('*')
            if p.is_file() and (p.name.endswith('.time.txt') or p.name == 'GNU_time.txt')}
    assert len(rows) == len(logs) == 584
    assert {r['log_path'] for r in rows} == logs
    assert len({r['log_path'] for r in rows}) == len(rows)
    for row in rows:
        assert (ROOT / row['record_path']).is_file(), row
        row['log_sha256'] = hashlib.sha256((ROOT / row['log_path']).read_bytes()).hexdigest()
    with (ROOT / 'results/TIMING_LOG_INDEX.tsv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], delimiter='\t')
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: row['log_path']))
    print('Indexed all 584 logs; no original log changed.')


if __name__ == '__main__':
    main()
