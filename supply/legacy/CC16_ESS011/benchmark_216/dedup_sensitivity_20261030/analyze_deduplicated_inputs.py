#!/usr/bin/env python3
"""Sensitivity analysis from saved timing rows and existing input SHA256 records.

No RDS is loaded and no inference, timing, or benchmark is performed.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import random
import statistics as st

import numpy as np
import pandas as pd


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path('/home/dzf/cellchat_acceleration/revision_2026_10')
    records = Path('/data/dzf/cellchatess_revision_20261026/results/common_external_timer_20261029/full')
    parser.add_argument('--timings', type=Path, default=records/'common_external_timer_216.tsv')
    parser.add_argument('--original-summary', type=Path, default=records/'common_external_timer_bootstrap_summary.tsv')
    parser.add_argument('--original-resamples', type=Path, default=records/'common_external_timer_bootstrap_resamples.tsv')
    parser.add_argument('--seed-audit', type=Path, default=root/'supplementary/tables/S2b_dataset_seed_audit.tsv')
    parser.add_argument('--input-hash-ledger', type=Path, default=root/'submission_materials_20261030/local_omitted_files.tsv')
    parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    timings = pd.read_csv(args.timings, sep='\t')
    audit = pd.read_csv(args.seed_audit, sep='\t')
    ledger = pd.read_csv(args.input_hash_ledger, sep='\t')
    hashes = ledger[ledger.logical_path.str.startswith('data/inputs/')].copy()
    hashes['experiment_id'] = hashes.logical_path.str.split('/').str[-2]
    hashes = hashes.rename(columns={'sha256': 'input_rds_sha256', 'local_path': 'input_rds_path'})
    joined = timings.merge(
        hashes[['experiment_id', 'input_rds_path', 'input_rds_sha256', 'hash_source']],
        on='experiment_id', validate='one_to_one'
    ).merge(
        audit[['dataset_id', 'target_cells', 'dataset_full_cells', 'three_seed_specific_downsamples']],
        on=['dataset_id', 'target_cells'], validate='many_to_one'
    )
    assert len(joined) == 216 and joined.dataset_id.nunique() == 12
    assert joined.CellChat_version.eq('1.6.1').all()
    assert joined.CellChatESS_version.eq('0.1.1').all()
    assert joined.timing_scope.eq('common_external_inference_call').all()
    assert np.allclose(joined.paired_speedup,
                       joined.original_inference_call_sec/joined.ess_inference_call_sec)
    assert joined.input_rds_sha256.str.fullmatch('[0-9a-f]{64}').all()
    joined['full_input'] = ~joined.three_seed_specific_downsamples
    full = joined[joined.full_input]
    assert len(full) == 99
    assert len(full[['dataset_id', 'target_cells']].drop_duplicates()) == 33
    assert full.groupby('dataset_id').input_rds_sha256.nunique().eq(1).all()
    assert joined.loc[~joined.full_input, 'input_rds_sha256'].nunique() == 117
    assert joined.input_rds_sha256.nunique() == 129

    dedup_rows = []
    for (dataset, digest), group in joined.groupby(['dataset_id', 'input_rds_sha256'], sort=True):
        assert group.actual_cells.nunique() == group.groups.nunique() == group.n_LR.nunique() == 1
        assert group.full_input.nunique() == 1
        assert len(group) == 1 or group.full_input.all()
        dedup_rows.append(dict(
            dataset_id=dataset, input_rds_sha256=digest,
            actual_cells=int(group.actual_cells.iloc[0]), K=int(group.groups.iloc[0]),
            n_LR=int(group.n_LR.iloc[0]), full_input=bool(group.full_input.iloc[0]),
            n_timing_records=len(group),
            target_scales=','.join(sorted(set(group.target_cells.astype(str)))),
            input_speedup=st.median(group.paired_speedup.tolist()),
            source_experiment_ids=';'.join(group.experiment_id),
        ))
    dedup = pd.DataFrame(dedup_rows)
    assert len(dedup) == 129 and dedup.n_timing_records.sum() == 216
    assert dedup.full_input.sum() == 12

    names = sorted(joined.dataset_id.unique())
    original_by = {ds: joined.loc[joined.dataset_id == ds, 'paired_speedup'].tolist() for ds in names}
    dedup_by = {ds: dedup.loc[dedup.dataset_id == ds, 'input_speedup'].tolist() for ds in names}
    rng = random.Random(20261026)
    bootstrap_rows = []
    for i in range(10000):
        sampled = [rng.choice(names) for _ in names]
        original_sample = [value for ds in sampled for value in original_by[ds]]
        dedup_sample = [value for ds in sampled for value in dedup_by[ds]]
        bootstrap_rows.append(dict(
            resample=i+1, sampled_datasets=';'.join(sampled),
            original_median=st.median(original_sample),
            deduplicated_median=st.median(dedup_sample),
            n_sampled_deduplicated_inputs=len(dedup_sample)
        ))
    bootstrap = pd.DataFrame(bootstrap_rows)
    reference_resamples = pd.read_csv(args.original_resamples, sep='\t')
    assert np.allclose(bootstrap.original_median, reference_resamples.paired_ratio_median,
                       rtol=1e-14, atol=1e-14)

    summaries = []
    for label, values, resample_col in [
        ('original_216_paired_timing_records', joined.paired_speedup.tolist(), 'original_median'),
        ('one_median_ratio_per_distinct_input_129', dedup.input_speedup.tolist(), 'deduplicated_median')
    ]:
        summaries.append(dict(
            analysis=label, n_records_or_inputs=len(values), n_datasets=len(names),
            median_speedup=st.median(values),
            ci95_lower=float(np.quantile(bootstrap[resample_col], .025, method='linear')),
            ci95_upper=float(np.quantile(bootstrap[resample_col], .975, method='linear')),
            bootstrap_replicates=10000, bootstrap_seed=20261026,
            bootstrap_unit='12 datasets sampled with replacement; retain all records/inputs in each selected dataset',
            ci_method='percentile 2.5/97.5; linear interpolation',
        ))
    summary = pd.DataFrame(summaries)
    ref_summary = pd.read_csv(args.original_summary, sep='\t').iloc[0]
    assert np.allclose(summary.loc[0, ['median_speedup','ci95_lower','ci95_upper']].to_numpy(dtype=float),
                       ref_summary[['estimate','ci95_lower','ci95_upper']].to_numpy(dtype=float))
    counts = joined.groupby('dataset_id').agg(
        original_timing_records=('experiment_id', 'size'),
        full_input_timing_records=('full_input', 'sum'),
    ).join(dedup.groupby('dataset_id').agg(
        distinct_inputs=('input_rds_sha256', 'size'),
        full_inputs=('full_input', 'sum'),
    )).reset_index()
    for frame, filename in [
        (joined, 'run_to_distinct_input_216.tsv'),
        (dedup, 'deduplicated_inputs_129.tsv'),
        (counts, 'dataset_input_counts_12.tsv'),
        (bootstrap, 'dataset_cluster_bootstrap_resamples_10000.tsv'),
        (summary, 'dedup_sensitivity_summary.tsv'),
    ]:
        frame.to_csv(args.out/filename, sep='\t', index=False, float_format='%.17g')

    delta = summary.median_speedup.iloc[1]-summary.median_speedup.iloc[0]
    percent = delta/summary.median_speedup.iloc[0]*100
    a, b = summary.iloc[0], summary.iloc[1]
    report = f'''# 统一计时 216 组的输入去重敏感性分析

只读取已有 TSV 和已有输入 SHA256 索引；不读取大型 RDS、不运行推断、不重测计时。

## 去重定义

先计算每条已保存配对记录的 Original/ESS 推断调用耗时比，再以
`dataset_id + input_rds_sha256` 标识同一输入，将其所有计时比取中位数。
没有用两种方法各自的耗时中位数相除。

33 个数据集–目标规模组合包含 99 条完整输入重复计时；跨规模核对 SHA256 后，
实际只有 12 份完整输入（每个数据集一份）。其他 117 条对应 117 份不同的下采样输入。
因此去重后为 **216 − 99 + 12 = 129 份输入**。
只在每个数据集–规模组合内压缩会得到 150 条，但仍重复保留跨规模相同的完整输入，
不能称为每个不同输入一条。

## 结果

| 分析 | 记录/输入数 | 中位加速比 | 数据集聚类 bootstrap 95% CI |
|:-----|------------:|-----------:|:--------------------------------|
| 原分析（配对计时记录） | 216 | {a.median_speedup:.2f}× | {a.ci95_lower:.2f}–{a.ci95_upper:.2f}× |
| 每份不同输入只保留一个中位比值 | 129 | {b.median_speedup:.2f}× | {b.ci95_lower:.2f}–{b.ci95_upper:.2f}× |

中位加速比变化 {delta:+.6f}×（相对变化 {percent:+.3f}%）。
这一变化反映不同输入的权重重新分配；两套估计都由已有统一边界计时记录计算。
此敏感性分析不消除固定执行顺序、并发进程竞争或输入之间的生物学依赖。

## Bootstrap 方法

与原分析相同：12 个数据集有放回抽取 12 次，抽中一个数据集即保留其全部去重输入，
重复 10,000 次；Python `random.Random(20261026)`；统计量为所抽输入比值的总体中位数；
区间为 2.5%/97.5% percentile，使用线性插值。去重后各数据集输入数不等，因此每次
bootstrap 的合并输入数可以变化；不在数据集内额外重采样，也没有改成数据集中位数的中位数。
本脚本以完全相同的数据集抽样序列重算原 216 条分析，逐次核对 10,000 个原 bootstrap
统计量，并复现原始区间。

## 文件与字段

- `run_to_distinct_input_216.tsv`：全部原记录与 SHA256/完整输入标记的对应关系。
- `deduplicated_inputs_129.tsv`：`input_speedup` 是同一输入的 `paired_speedup` 中位数。
- `dataset_input_counts_12.tsv`：每个数据集去重前后的计数。
- `dataset_cluster_bootstrap_resamples_10000.tsv`：全部重采样结果和抽中的数据集 ID。
- `dedup_sensitivity_summary.tsv`：精确估计及区间。
- `source_checksums.tsv`：输入 TSV/脚本的绝对路径、大小及 SHA256。

计时环境：CellChat 1.6.1 / CellChatESS 0.1.1；计时字段为
`original_inference_call_sec` 和 `ess_inference_call_sec`，从紧邻函数调用前到返回后计时，
包含函数内部转换及初始化，排除启动、加载、I/O 和下游分析。
输入 SHA256 取自已有 `local_omitted_files.tsv` 中 `data/inputs/` 行的 `sha256` 列；
完整输入标记取自 `S2b_dataset_seed_audit.tsv` 的 `three_seed_specific_downsamples`。

## 可用于回复信的文字

为检验对同一完整输入的重复计时是否影响总体结果，我们按输入文件 SHA256 合并重复记录，
每份输入保留一个配对加速比中位数。99 条完整输入计时跨 33 个数据集–规模组合，实际对应
12 份完整输入；与 117 份不同的下采样输入合并后，共保留 129 份不同输入。
去重后的中位加速比为 {b.median_speedup:.2f}×（95% 数据集聚类 bootstrap CI：
{b.ci95_lower:.2f}–{b.ci95_upper:.2f}×；10,000 次重采样），相比原分析的
70.82×（52.35–93.46×），变化为 {percent:+.1f}%。
'''
    (args.out/'dedup_sensitivity_report_zh.md').write_text(report)
    sources = [args.timings, args.original_summary, args.original_resamples,
               args.seed_audit, args.input_hash_ledger, Path(__file__).resolve()]
    pd.DataFrame([dict(path=str(p.resolve()), size_bytes=p.stat().st_size, sha256=sha256(p))
                  for p in sources]).to_csv(args.out/'source_checksums.tsv', sep='\t', index=False)
    (args.out/'analysis_environment.json').write_text(json.dumps(dict(
        python=platform.python_version(), pandas=pd.__version__, numpy=np.__version__,
        no_inference_or_timing_rerun=True, original_bootstrap_all_10000_resamples_verified=True,
    ), indent=2)+'\n')
    print(summary[['analysis','n_records_or_inputs','median_speedup','ci95_lower','ci95_upper']].to_string(index=False))
    print(f'Change: {delta:+.9f}x ({percent:+.6f}%)')


if __name__ == '__main__':
    main()
