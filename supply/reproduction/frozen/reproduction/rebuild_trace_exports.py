#!/usr/bin/env python3
"""Rejoin direct trace files without R's automatic hex-to-numeric conversion.

This performs no inference and leaves all executed-run records unchanged.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--delivery', type=Path, required=True)
a = p.parse_args()
root = a.delivery.resolve()

def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def key(row):
    return tuple(int(row[k]) for k in ('lr_index_1based', 'pair_index_0based', 'bootstrap'))

def yes(value):
    return str(value).lower() in ('true', '1')

report = []
for version in ('cc16', 'cc22'):
    folder = root / 'results/new/xenium_target_trace' / version / 'output'
    paths = [folder / 'original_direct_trace.tsv', folder / 'ess_direct_trace.tsv']
    original, ess = (read(path) for path in paths)
    omap, emap = ({key(r): r for r in rows} for rows in (original, ess))
    assert len(original) == len(omap) == len(ess) == len(emap) == 500
    assert omap.keys() == emap.keys()
    out = []
    for k in sorted(omap):
        row = dict(zip(('lr_index_1based', 'pair_index_0based', 'bootstrap'), k))
        decisions = []
        for label, data in (('original', omap[k]), ('ESS', emap[k])):
            obs, boot = float(data['P_obs']), float(data['P_boot'])
            assert obs == float.fromhex(data['P_obs_hex'])
            assert boot == float.fromhex(data['P_boot_hex'])
            assert float(data['delta']) == boot - obs
            decision = boot > obs
            assert decision == yes(data['reject'])
            decisions.append(decision)
            for field in ('P_obs', 'P_boot', 'delta'):
                row[label + '_' + field] = data[field]
            row[label + '_P_obs_hex'] = data['P_obs_hex']
            row[label + '_P_boot_hex'] = data['P_boot_hex']
            row[label + '_delta_hex'] = float(data['delta']).hex()
            row[label + '_reject'] = decision
        row['decision_equal'] = decisions[0] == decisions[1]
        out.append(row)
    summary = read(folder / 'target_summary.tsv')
    for target in summary:
        chosen = [r for r in out if r['lr_index_1based'] == int(target['lr1']) and r['pair_index_0based'] == int(target['pair0'])]
        assert len(chosen) == 100
        for label in ('original', 'ESS'):
            count = sum(r[label + '_reject'] for r in chosen)
            assert count == int(target[label + '_count'])
            assert count / 100 == float(target[label + '_pvalue'])
    dest = root / 'supplementary' / ('S_G2_direct_decisions_' + version + '.tsv')
    with dest.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter='\t')
        w.writeheader()
        w.writerows(out)
    report.append(dict(environment=version, comparisons=len(out), targets=len(summary),
                       different_decisions=sum(not r['decision_equal'] for r in out),
                       source_files=[str(path.relative_to(root)) for path in paths],
                       source_sha256=[hashlib.sha256(path.read_bytes()).hexdigest() for path in paths],
                       output=str(dest.relative_to(root)), hex_decimal_roundtrip='PASS',
                       strict_comparison_recheck='PASS', final_counts_and_pvalues='PASS'))
(root / 'audit/trace_export_validation.json').write_text(json.dumps(report, indent=2) + '\n')
(root / 'audit/trace_export_note_CN.md').write_text('''# 十六进制追踪导出的整理说明

实际 C++ 直接追踪文件 ess_direct_trace.tsv 保留了正确的十六进制浮点字符串。
运行时 run_target_trace_CC22.R 第80行的 read.delim(trace_path) 自动将这些字符串转成数值；
因此旧 comparison_all_decisions.tsv 中 ESS_Pobs_hex / ESS_Pboot_hex 两列误导出成十进制。
这只是汇合表的导出类型问题，未影响实际追踪、严格比较、计数或 p 值。

本包不覆盖原运行文件。supplementary/S_G2_direct_decisions_cc16.tsv 和 cc22.tsv
用 reproduction/rebuild_trace_exports.py 从两侧原始直接追踪文件重新汇合，保留原始17位十进制和十六进制。
脚本验证每个浮点数十进制/十六进制可逆、delta、严格 > 判定及最终计数/p值，共500条/环境。
没有重新运行推断，也没有改变概率或 p 值公式。

未来执行入口的最小整理修复见 trace_export_colclasses.patch；该补丁尚未应用到本次实际运行脚本，
从而保留运行时源码身份。不能把重新汇合表说成一批新的实验。
''')
(root / 'audit/trace_export_colclasses.patch').write_text('''--- a/reproduction/run_target_trace_CC22.R
+++ b/reproduction/run_target_trace_CC22.R
@@ -80 +80 @@
-o<-do.call(rbind,original_records);e<-read.delim(trace_path)
+o<-do.call(rbind,original_records);e<-read.delim(trace_path,colClasses=c(P_obs_hex="character",P_boot_hex="character"))
''')
print(json.dumps(report, indent=2))
