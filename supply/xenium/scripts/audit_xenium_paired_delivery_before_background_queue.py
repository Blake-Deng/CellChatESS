#!/usr/bin/env python3
"""Index terminal strict-pair evidence; never run inference or relabel a run."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def read_json(path):
    return json.loads(path.read_text()) if path.is_file() else {}

def write_tsv(path, fields, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    root = args.output.resolve()
    pipeline = read_json(root / 'pipeline_status.json')
    if pipeline.get('status') == 'RUNNING':
        raise SystemExit('Pipeline still RUNNING: refusing a terminal delivery audit')
    claims = []
    methods = {}
    for name in ('original', 'ess'):
        status = read_json(root / name / 'run_status.json')
        meta = read_json(root / name / 'run_metadata.json')
        complete = status.get('status') == 'COMPLETED' and bool(meta)
        path = root / name / 'result.rds'
        hash_ok = complete and path.is_file() and sha(path) == meta['outputs']['output_rds_sha256']
        methods[name] = {'status': status.get('status', 'NOT_RUN'),
                         'complete': complete, 'result_hash_matches': hash_ok,
                         'evidence': str(Path(name) / 'run_status.json')}
        claims.append({'claim': name + ' full 627-LR output',
                       'status': 'PASS' if hash_ok else methods[name]['status'],
                       'value_or_scope': 'This frozen Xenium input only',
                       'evidence': methods[name]['evidence']})
    summary_path = root / 'comparison/full627_comparison_summary.tsv'
    comparison = read_json(root / 'comparison_status.json')
    summaries = {}
    if summary_path.is_file():
        with summary_path.open() as handle:
            summaries = {r['check']: r['value'] for r in csv.DictReader(handle, delimiter='\t')}
    full_pair = all(x['result_hash_matches'] for x in methods.values())
    full_comparison = full_pair and comparison.get('status') == 'COMPLETED' and bool(summaries)
    for key in ('maximum_absolute_probability_difference', 'probability_different_elements',
                'pvalue_different_elements', 'aggregate_counts_identical',
                'original_aggregate_count_recomputed_matches_saved',
                'ESS_aggregate_count_recomputed_matches_saved'):
        claims.append({'claim': key, 'status': 'MEASURED' if full_comparison else 'NOT_CONFIRMED',
                       'value_or_scope': summaries.get(key, '无法确认') if full_comparison else '无法确认',
                       'evidence': 'comparison/full627_comparison_summary.tsv'})
    write_tsv(root / 'terminal_claim_evidence_map.tsv',
              ['claim', 'status', 'value_or_scope', 'evidence'], claims)
    audit = {'pipeline': pipeline, 'methods': methods,
             'both_complete_and_output_hashes_verified': full_pair,
             'full_comparison_completed': full_comparison,
             'scope': 'Only this frozen input, source builds and actual conditions; no universal identity claim',
             'large_file_policy': 'Large original outputs remain in place; hashes and paths are indexed, not omitted from verification'}
    (root / 'terminal_delivery_integrity.json').write_text(json.dumps(audit, indent=2, ensure_ascii=False) + '\n')
    lines = ['# 严格配对验证交付索引', '',
             '本索引仅在队列进入终止状态后生成，不运行推断，不修改原始日志或软件源码。', '',
             '| 项目 | 状态 | 证据 |', '|---|---|---|']
    for name, result in methods.items():
        lines.append(f"| {name} 全量输出 | {result['status']}；输出哈希核对 {result['result_hash_matches']} | {result['evidence']} |")
    lines += ['', '完整配对和逐元素比较已完成：' + ('是' if full_comparison else '否') + '。']
    if full_comparison:
        lines += ['', '完整中文结果：`Xenium_627LR_paired_report_CN.md`。',
                  '并列表 R2：`R2_xenium_memory_time_comparison.tsv`；同口径比值：`paired_same_scope_resource_ratios.tsv`。',
                  '全部数值、非有限值和原五目标：`comparison/`；差异直接追踪：`difference_diagnostics/`。',
                  'RSS 原始采样：各方法 `rss_samples.tsv`；可编辑资源图：`xenium_paired_RSS.*`。']
    else:
        lines += ['', '无法确认全量数值一致性；不得计算或引用完整配对性能比。',
                  '失败/中止阶段以各方法 run_status.json、run.log、GNU time 和 stage_events.tsv 为准。',
                  'smoke、局部 TRACE、旧 PROFILE 均不替代未完成的 627-LR 全量结果。']
    lines += ['', '源码与输入身份：`source_and_install_identity.json`、各方法 `run_identity_before.json`；',
              '源码/归档：`sources/`；环境锁：`environment/`；脚本及五目标配置：`scripts/`。',
              '对象大小、纯载荷和 RSS 分开报告；每个稠密输出仍为 O(nLR K²)。',
              '未对干净计时构建的内部临时对象作插桩。源码推导的载荷上界不等于实测同时存活工作区。',
              '终止后的静态文件 SHA256：`terminal_SHA256_manifest.tsv`。动态日志和状态文件另列，不冒称已冻结。', '']
    (root / 'DELIVERY_INDEX_CN.md').write_text('\n'.join(lines))
    excluded = {'terminal_SHA256_manifest.tsv', 'SHA256_manifest.tsv', 'delivery_audit.log',
                'finalization.log', 'pipeline.log', 'pipeline_status.json'}
    entries = []
    large = []
    for path in sorted(root.rglob('*')):
        if path.is_file() and path.name not in excluded:
            row = {'relative_path': str(path.relative_to(root)),
                   'size_bytes': path.stat().st_size, 'sha256': sha(path)}
            entries.append(row)
            if row['size_bytes'] > 50 * 1024 ** 2:
                large.append({**row, 'absolute_path': str(path)})
    write_tsv(root / 'large_outputs_manifest.tsv',
              ['relative_path', 'size_bytes', 'sha256', 'absolute_path'], large)
    entries.append({'relative_path': 'large_outputs_manifest.tsv',
                    'size_bytes': (root / 'large_outputs_manifest.tsv').stat().st_size,
                    'sha256': sha(root / 'large_outputs_manifest.tsv')})
    # Avoid a stale checksum if this independent audit is deliberately repeated.
    entries = [r for r in entries if r['relative_path'] != 'large_outputs_manifest.tsv'] + [entries[-1]]
    write_tsv(root / 'terminal_SHA256_manifest.tsv', ['relative_path', 'size_bytes', 'sha256'], entries)
    print(json.dumps(audit, indent=2, ensure_ascii=False))

if __name__ == '__main__':
    main()
