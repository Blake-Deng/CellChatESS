#!/usr/bin/env python3
"""Refresh distributable-file hashes after reviewing generated changes."""
from pathlib import Path
import csv
import hashlib

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / 'MANIFEST.tsv'

if __name__ == '__main__':
    files = sorted(p for p in ROOT.rglob('*') if p.is_file() and p != TARGET
                   and '.git' not in p.parts and '__pycache__' not in p.parts)
    with TARGET.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=['path', 'size_bytes', 'sha256'], delimiter='\t')
        writer.writeheader()
        for path in files:
            writer.writerow({'path': str(path.relative_to(ROOT)),
                             'size_bytes': path.stat().st_size,
                             'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    print(f'Manifest refreshed: {len(files)} files.')
