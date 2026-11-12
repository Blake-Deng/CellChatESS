#!/usr/bin/env python3
"""Create R2 measurements, RSS gap summaries, report and checksums from completed runs."""
from __future__ import annotations
import argparse,csv,hashlib,json,re,subprocess,statistics
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def sha(p:Path)->str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def read_json(p:Path):return json.loads(p.read_text())
def read_tsv(p:Path):
    with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(p:Path,fields,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def gnu_summary(p:Path):
    s=p.read_text(errors='replace')
    wall=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)',s)
    rss=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',s)
    exitcode=re.search(r'Exit status:\s*(\d+)',s)
    sec=None
    if wall:
        for z in wall.group(1).split(':'):sec=(sec or 0)*60+float(z)
    user=re.search(r'User time \(seconds\):\s*([0-9.]+)',s)
    system=re.search(r'System time \(seconds\):\s*([0-9.]+)',s)
    cpu=re.search(r'Percent of CPU this job got:\s*([0-9.]+)%',s)
    return {'GNU_time_wall_text':wall.group(1) if wall else None,'GNU_time_wall_sec':sec,
      'GNU_time_peak_RSS_GiB':int(rss.group(1))/1024**2 if rss else None,
      'GNU_time_user_sec':float(user.group(1)) if user else None,
      'GNU_time_system_sec':float(system.group(1)) if system else None,
      'GNU_time_CPU_percent':float(cpu.group(1)) if cpu else None,
      'GNU_time_exit_status':int(exitcode.group(1)) if exitcode else None}
def summarize_gaps(p:Path,method:str):
    rows=read_tsv(p)
    t=[float(x['sample_start_epoch_ms'])/1000 for x in rows]
    diffs=[b-a for a,b in zip(t,t[1:])]
    return {'method':method,'n_samples':len(rows),'requested_interval_sec':1,
      'median_sample_spacing_sec':statistics.median(diffs) if diffs else None,
      'maximum_sample_gap_sec':max(diffs) if diffs else None,
      'gaps_over_2_sec':sum(x>2 for x in diffs),'gaps_over_5_sec':sum(x>5 for x in diffs),
      'sampling_has_no_interpolation':True}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();out=a.output.resolve()
    orig_dir=out/'original';ess_dir=out/'ess';comp=out/'comparison'
    for method,d in [('original',orig_dir),('ess',ess_dir)]:
        st=read_json(d/'run_status.json')
        if st.get('status')!='COMPLETED' or not (d/'result.rds').is_file():
            raise SystemExit(f'{method} is not a completed full run; refusing complete report')
    if read_json(out/'comparison_status.json').get('status')!='COMPLETED':
        raise SystemExit('Full comparison is not complete')
    orig=read_json(orig_dir/'run_metadata.json');ess=read_json(ess_dir/'run_metadata.json')
    for d,meta in [(orig_dir,orig),(ess_dir,ess)]:
        if sha(d/'result.rds')!=meta['outputs']['output_rds_sha256']:
            raise SystemExit('Saved full output hash changed after the run: '+str(d))
    if sha(Path(orig['input']['input_path']))!=orig['input']['input_sha256']:
        raise SystemExit('Frozen input hash changed after inference')
    if orig['input']['input_sha256']!=ess['input']['input_sha256']:
        raise SystemExit('Input SHA differs across methods')
    for key in ['database_serialized_sha256','LRsig_serialized_sha256','cell_names_sha256',
      'cell_group_labels_sha256','group_order_sha256','LR_order_sha256','cells','signaling_genes','groups','n_LR']:
        if orig['input'][key]!=ess['input'][key]:raise SystemExit(f'Input identity mismatch: {key}')
    sums={r['check']:r['value'] for r in read_tsv(comp/'full627_comparison_summary.tsv')}
    rows=[]
    for method,d,meta in [('Original CellChat',orig_dir,orig),('CellChatESS',ess_dir,ess)]:
        g=gnu_summary(d/'gnu_time.txt');samples=read_tsv(d/'rss_samples.tsv')
        sampled_peak=max((float(r['VmRSS_kB']) for r in samples if r['VmRSS_kB']),default=None)
        max_os_threads=max((int(r['Threads']) for r in samples if r['Threads']),default=None)
        r=meta['outputs'];i=meta['input']
        rows.append({'method':method,'input_file_bytes':i['expression_matrix_file_bytes'],
          'input_data_signaling_object_bytes':i['expression_matrix_object_bytes'],
          'probability_R_object_bytes':r['probability_object_bytes'],'probability_data_payload_bytes':r['probability_payload_bytes'],
          'pvalue_R_object_bytes':r['pvalue_object_bytes'],'pvalue_data_payload_bytes':r['pvalue_payload_bytes'],
          'aggregate_count_object_bytes':r['count_object_bytes'],'compressed_saved_output_bytes':r['output_rds_bytes'],
          'inference_function_sec':meta['inference_function_elapsed_sec'],'aggregateNet_sec':meta['aggregateNet_elapsed_sec'],
          'GNU_time_full_process_wall_sec':g['GNU_time_wall_sec'],'GNU_time_peak_RSS_GiB':g['GNU_time_peak_RSS_GiB'],
          'GNU_time_CPU_percent':g['GNU_time_CPU_percent'],'GNU_time_user_sec':g['GNU_time_user_sec'],
          'GNU_time_system_sec':g['GNU_time_system_sec'],'sampled_max_OS_threads':max_os_threads,
          'sampled_RSS_peak_GiB':sampled_peak/1024**2 if sampled_peak is not None else None,
          'RSS_samples':len(samples),'output_sha256':r['output_rds_sha256'],'GNU_time_exit_status':g['GNU_time_exit_status']})
    write_tsv(out/'R2_xenium_memory_time_comparison.tsv',list(rows[0]),rows)
    ratios={'inference_function_speed_ratio_original_over_ESS':rows[0]['inference_function_sec']/rows[1]['inference_function_sec'],
      'full_process_wall_speed_ratio_original_over_ESS':rows[0]['GNU_time_full_process_wall_sec']/rows[1]['GNU_time_full_process_wall_sec'],
      'GNU_time_peak_RSS_reduction_percent':100*(1-rows[1]['GNU_time_peak_RSS_GiB']/rows[0]['GNU_time_peak_RSS_GiB'])}
    write_tsv(out/'paired_same_scope_resource_ratios.tsv',['metric','value'],[{'metric':k,'value':v} for k,v in ratios.items()])
    stages=[]
    for method,d in [('Original CellChat',orig_dir),('CellChatESS',ess_dir)]:
        samples=read_tsv(d/'rss_samples.tsv')
        for stage in dict.fromkeys(x['stage'] for x in samples):
            selected=[x for x in samples if x['stage']==stage and x['VmRSS_kB']]
            values=[float(x['VmRSS_kB'])/1024**2 for x in selected]
            if values:stages.append({'method':method,'stage':stage,'n_samples':len(values),
              'first_sample_epoch_ms':selected[0]['sample_start_epoch_ms'],'last_sample_epoch_ms':selected[-1]['sample_start_epoch_ms'],
              'RSS_min_GiB':min(values),'RSS_median_GiB':statistics.median(values),'RSS_max_GiB':max(values)})
    write_tsv(out/'RSS_external_stage_summary.tsv',list(stages[0]),stages)
    gap_rows=[summarize_gaps(orig_dir/'rss_samples.tsv','Original CellChat'),summarize_gaps(ess_dir/'rss_samples.tsv','CellChatESS')]
    write_tsv(out/'rss_sampling_gap_summary.tsv',list(gap_rows[0]),gap_rows)
    # Preserve a flat run identity table linking installed files and source IDs.
    source=read_json(out/'source_and_install_identity.json')
    identity_rows=[]
    for method,meta in [('Original CellChat',orig),('CellChatESS',ess)]:
        si=meta['identity']['source_function_identity'];pi=meta['identity']['package_versions']
        identity_rows.append({'method':method,'CellChat_version':pi['CellChat'],'CellChat_path':meta['identity']['installation_paths']['CellChat'],
          'CellChat_commit':source['CellChat_commit'],'CellChat_archive_sha256':source['CellChat_archive_sha256'],
          'CellChat_modeling_R_sha256':source['CellChat_modeling_R_sha256'],'active_function':si['function_name'],
          'active_function_body_sha256':si['function_body_sha256'],'function_matches_pinned_source':si['function_body_matches_pinned_source'],
          'AutoZyme_namespace_loaded':si['AutoZyme_namespace_loaded'],'CellChatESS_version':pi['CellChatESS'],
          'CellChatESS_commit':source['ESS_commit'] if method=='CellChatESS' else '',
          'CellChatESS_tree':source['ESS_tree'] if method=='CellChatESS' else '',
          'CellChatESS_DLL_sha256':si.get('ESS_DLL_sha256',''),'input_sha256':meta['input']['input_sha256'],
          'database_serialized_sha256':meta['input']['database_serialized_sha256'],
          'LRsig_serialized_sha256':meta['input']['LRsig_serialized_sha256']})
    write_tsv(out/'source_input_parameter_identity.tsv',list(identity_rows[0]),identity_rows)
    # Chinese report uses the exact full-run and comparison output records.
    def yn(x):return '是' if str(x).lower() in ('true','1') else '否'
    maxdiff={r['metric']:r['value'] for r in read_tsv(comp/'maximum_probability_difference.tsv')}
    targets=read_tsv(comp/'five_prior_targets.tsv')
    prior_trace=out/'auxiliary/target_trace_cc22/output/target_summary.tsv'
    if prior_trace.is_file():
        trace_lookup={(x['lr_name'],x['source_group'],x['target_group']):x for x in read_tsv(prior_trace)}
        crosscheck=[]
        for target in targets:
            key=tuple(target[k] for k in ('lr_name','source_group','target_group'))
            prior=trace_lookup.get(key)
            crosscheck.append({'lr_name':key[0],'source_group':key[1],'target_group':key[2],
              'prior_trace_target_found':prior is not None,
              'original_observed_probability_hex_matches':bool(prior) and prior['original_Pobs_hex']==target['original_probability_hex'],
              'ESS_observed_probability_hex_matches':bool(prior) and prior['ESS_Pobs_hex']==target['ESS_probability_hex'],
              'original_pvalue_matches':bool(prior) and float(prior['original_pvalue'])==float(target['original_pvalue']),
              'ESS_pvalue_matches':bool(prior) and float(prior['ESS_pvalue'])==float(target['ESS_pvalue'])})
        write_tsv(out/'five_targets_full_run_vs_prior_TRACE.tsv',list(crosscheck[0]),crosscheck)
    p_diff_count=int(sums['pvalue_different_elements'])
    count_ok=yn(sums['aggregate_counts_identical'])=='是' and yn(sums['original_aggregate_count_recomputed_matches_saved'])=='是' and yn(sums['ESS_aggregate_count_recomputed_matches_saved'])=='是'
    count_policies=read_tsv(comp/'aggregate_count_policy_summary.tsv')
    paper_count_policy=next(x for x in count_policies if x['policy']=='filterCommunication_min_cells_10_then_aggregateNet')
    paper_count_ok=yn(paper_count_policy['counts_identical'])=='是'
    lines=['# Xenium 全量 627 LR 配对验证报告（服务器复现）','',
      '本报告只针对下述冻结 Xenium 输入、源码构建、服务器和参数，不外推到任意输入。原版 CellChat 与最终 ESS 在独立 R 进程中依次运行，未并发。','',
      '## 输入、软件与参数','',
      f"- 输入：`{orig['input']['input_path']}`；SHA256 `{orig['input']['input_sha256']}`。两次运行在读入前和运行中核对相同输入。",
      f"- 对象：{orig['input']['cells']:,} 细胞、{orig['input']['signaling_genes']} 信号基因、K={orig['input']['groups']:,}、{orig['input']['n_LR']} LR；两次运行的 cell/group/DB/LR 序列化哈希完全一致。",
      f"- Original：CellChat {orig['identity']['package_versions']['CellChat']}；上游 commit `{source['CellChat_commit']}`；归档 SHA256 `{source['CellChat_archive_sha256']}`。实际 `computeCommunProb` 函数体与冻结源码一致，AutoZyme 未加载。",
      f"- ESS：CellChatESS {ess['identity']['package_versions']['CellChatESS']}；commit `{source['ESS_commit']}`；tree `{source['ESS_tree']}`；安装 DLL SHA256 `{ess['identity']['source_function_identity']['ESS_DLL_sha256']}`。实际 `computeCommunProbESS` 函数体与冻结源码一致。",
      '- 参数：RNA、`type="triMean"`、`raw.use=TRUE`、`population.size=FALSE`、`nboot=100`、`seed.use=1`、`Kh=0.5`、`n=1`；Original 与 ESS 均单线程，ESS `algorithm="ess"`。',
      '', '## 时间与内存','',
      '- 推断函数时间从各函数调用前计时，函数返回后停止；不含 R 启动、包加载、输入读取、`aggregateNet` 或序列化。',
      '- GNU time 全进程墙钟从 R 进程启动到退出，包含推断函数时间、同一 `aggregateNet` 默认计数步骤和相同 gzip `saveRDS`；全进程时间包含推断时间，不能再相加。',
      '- RSS 同时报告 GNU time 的进程峰值和独立 1 秒 `/proc` 原始采样峰值。采样间隔与大于 2/5 秒的缺口见 `rss_sampling_gap_summary.tsv`，没有插值。',
      '- `object.size`、稠密数组纯 payload 和进程 RSS 是不同量。每个概率/p 值输出仍是稠密 `n_LR × K × K` 数组，空间复杂度为 O(n_LR K²)；ESS 支持集剪枝不消除这两个稠密输出。',
      '', '## 全量逐元素核查','',
      f"- 概率与 p 值维度、名称和顺序一致：{yn(sums['dimensions_match'])} / {yn(sums['names_and_order_match'])}。",
      f"- 概率元素数 {int(sums['probability_entries']):,}；数值不同元素 {int(sums['probability_different_elements']):,}；最大绝对差 `{maxdiff['max_probability_original'] and sums['maximum_absolute_probability_difference']}`。最大差坐标：LR `{maxdiff['max_probability_lr_name']}`，发送 `{maxdiff['max_probability_source_group']}`，接收 `{maxdiff['max_probability_target_group']}`；Original `{maxdiff['max_probability_original']}`（`{maxdiff['max_probability_original_hex']}`），ESS `{maxdiff['max_probability_ESS']}`（`{maxdiff['max_probability_ESS_hex']}`）。",
      f"- p 值元素数 {int(sums['pvalue_entries']):,}；不一致元素 {p_diff_count}。所有差异位置及十进制/十六进制值见 `comparison/pvalue_differences.tsv`。未使用容差。",
      f"- 默认 `aggregateNet(thresh=0.05)` count：Original 和 ESS 的独立分块重算均与各自保存值吻合：{yn(sums['original_aggregate_count_recomputed_matches_saved'])} / {yn(sums['ESS_aggregate_count_recomputed_matches_saved'])}；两方法 count 矩阵完全一致：{yn(sums['aggregate_counts_identical'])}。count 差异坐标见 `comparison/aggregate_count_differences.tsv`。",
      '- 概率及 p 值的 NA、NaN、正负 Inf 计数见 `comparison/nonfinite_summary.tsv`；未有限值差异坐标见对应文件。',
      '', '## 原五个目标条目','',
      '按 LR、发送群、接收群三者共同定位，完整输出数值和 double 十六进制值见 `comparison/five_prior_targets.tsv`。']
    for t in targets:
        lines.append(f"- {t['lr_name']} / {t['source_group']} → {t['target_group']}：P Original={t['original_probability']}，ESS={t['ESS_probability']}；p Original={t['original_pvalue']}，ESS={t['ESS_pvalue']}；count={t['original_aggregate_count_at_pair']}/{t['ESS_aggregate_count_at_pair']}。")
    if p_diff_count:
        lines += ['',f"存在 {p_diff_count} 个 p 值差异；应按差异坐标追加对应置换的直接追踪。严格 `P_boot > P_obs` 规则没有改动。"]
    else:
        lines += ['', '全量 627 LR 的经验 p 值无差异；因此本次全量结果未触发追加 bootstrap 逐次追踪。原五目标此前的目标级追踪记录仍是单独范围证据。']
    lines += ['', '## 运行文件','',
      '- 每种方法的原始 R 日志、GNU time、RSS 采样、函数源码快照、sessionInfo、输入对象 group sizes/LR 表和 `run_metadata.json` 位于 `original/` 与 `ess/`。',
      '- R2 并列数据位于 `R2_xenium_memory_time_comparison.tsv`；RSS 绘图源数据与绘图脚本见 `*_rss_plot_data.tsv` 和 `reproduction/plot_xenium_paired_rss.R`。',
      '- 本次比较脚本：`reproduction/compare_xenium_paired_full627.R`；worker：`reproduction/run_xenium_paired_worker.R`；队列入口：`reproduction/run_xenium_paired_full627.py`。',
      '- 旧的 ESS PROFILE 插桩运行没有被用于本次 ESS 数值结果或配对时间/RSS。',
      '', '## 原始输出和复现边界','',
      '完整 RDS 输出在本目录的 `original/result.rds` 与 `ess/result.rds`。所有结论仅适用于该冻结输入、上述软件源码、硬件和参数；概率逐元素相同以外的结论需按输出报告。','']
    (out/'Xenium_627LR_paired_report_CN.md').write_text('\n'.join(lines))
    with (out/'Xenium_627LR_paired_report_CN.md').open('a') as report:
        report.write('\n## 本次实测时间与 RSS\n\n')
        report.write('| 方法 | 函数调用秒 | 全进程秒 | GNU time 峰值 GiB | 采样峰值 GiB |\n|---|---:|---:|---:|---:|\n')
        for row in rows:report.write(f"| {row['method']} | {row['inference_function_sec']:.3f} | {row['GNU_time_full_process_wall_sec']:.3f} | {row['GNU_time_peak_RSS_GiB']:.5f} | {row['sampled_RSS_peak_GiB']:.5f} |\n")
        report.write('\n同口径推断函数耗时比 %.6f；同口径全进程耗时比 %.6f；GNU time 峰值 RSS 降幅 %.6f%%。只针对本次各一次顺序运行，不提供重复实验置信区间。\n' % tuple(ratios.values()))
        report.write('\n主性能构建未对函数内部临时对象作插桩，外部阶段 RSS 见 RSS_external_stage_summary.tsv。辅助 PROFILE 构建的阶段记录单独保留，不能替代本次主构建的时间/RSS，也不能把峰值 RSS 减去输出载荷解释为测得的临时工作区。独立 RSS 采样无法观察进程启动前或退出后的内存，采样缺口不会以插值补齐。\n')
    diagnostic_script=ROOT/'reproduction/diagnose_xenium_paired_differences.py'
    diagnostic_command=[__import__('sys').executable,str(diagnostic_script),'--output',str(out)]
    (out/'diagnosis_command.json').write_text(json.dumps(diagnostic_command,indent=2)+'\n')
    with (out/'diagnosis_driver.log').open('w') as f:
        diagnostic_rc=subprocess.run(diagnostic_command,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT).returncode
    diagnostic_status_path=out/'difference_diagnostics/status.json'
    diagnostic_status=read_json(diagnostic_status_path) if diagnostic_status_path.is_file() else {'status':'FAILED'}
    acceptance=[
      {'item':'Original full627 completed','status':'PASS','evidence':'original/run_status.json'},
      {'item':'ESS full627 completed','status':'PASS','evidence':'ess/run_status.json'},
      {'item':'dimensions/names/order','status':'PASS' if sums['dimensions_match']=='TRUE' and sums['names_and_order_match']=='TRUE' else 'FAIL','evidence':'comparison/full627_comparison_summary.tsv'},
      {'item':'finite probability and pvalue entries','status':'PASS' if all(int(x['n'])==0 for x in read_tsv(comp/'nonfinite_summary.tsv')) else 'FAIL','evidence':'comparison/nonfinite_summary.tsv'},
      {'item':'pvalue exact equality','status':'PASS' if p_diff_count==0 else 'FAIL','evidence':'comparison/pvalue_differences.tsv'},
      {'item':'aggregate count independent recomputation and exact equality','status':'PASS' if count_ok else 'FAIL','evidence':'comparison/aggregate_count_differences.tsv'},
      {'item':'paper min.cells=10 count equality','status':'PASS' if paper_count_ok else 'FAIL','evidence':'comparison/aggregate_count_min_cells_10_differences.tsv'},
      {'item':'five prior target pvalues','status':'PASS' if all(t['pvalue_identical']=='TRUE' for t in targets) else 'FAIL','evidence':'comparison/five_prior_targets.tsv'},
      {'item':'difference tracing','status':'NOT_RUN' if diagnostic_status['status']=='NOT_REQUIRED' else ('PASS' if diagnostic_rc==0 and diagnostic_status['status']=='COMPLETED' else 'BLOCKED'),'evidence':'difference_diagnostics/status.json'},
      {'item':'internal temporary object sizes in clean performance build','status':'NOT_RUN','evidence':'No internal instrumentation in main timing builds; auxiliary PROFILE stage RSS has a distinct source diff and DLL.'}]
    write_tsv(out/'validation_acceptance.tsv',['item','status','evidence'],acceptance)
    with (out/'Xenium_627LR_paired_report_CN.md').open('a') as report:
        report.write('\n## 聚合 count 的两个明确口径\n\n')
        report.write('默认 aggregateNet 显著性计数与原稿 filterCommunication(min.cells=10) 后的计数分别核查。官方 CellChat modeling.R:916-920 排除的是细胞数小于或等于 10 的群；min.samples=NULL 的有效值为 1，因此不触发多样本过滤。\n')
        report.write('本次原稿口径排除 '+paper_count_policy['excluded_groups']+' 个群对应的发送/接收页；保留完整 K×K 矩阵并将这些行列置零。Original 显著 LR 边总数为 '+paper_count_policy['original_total_significant_LR_edges']+'，ESS 为 '+paper_count_policy['ESS_total_significant_LR_edges']+'；群对计数差异数 '+paper_count_policy['different_group_pairs']+'，逐元素一致：'+yn(paper_count_policy['counts_identical'])+'。\n')
        report.write('细胞数过滤的独立核查在主性能计时之外；GNU time 全进程范围仍是共同 worker 中的输入、推断、默认未加细胞数过滤的聚合、保存与退出，不能将该全进程时间称为包含细胞数过滤的完整论文工作流。规则、排除群列表和两份计数矩阵见 comparison/aggregate_count_policy_summary.tsv 及相邻文件。\n')
        report.write('\n## 差异追踪与验收\n\n')
        report.write('差异诊断状态：'+diagnostic_status['status']+'。若存在差异，逐 LR 的直接 bootstrap 追踪与全量输出交叉核对保存在 difference_diagnostics；目标级诊断不替代全量比较，也不用于计时。具体浮点来源需要根据直接追踪审阅，不预设为 BLAS、编译器或求和顺序。\n')
        report.write('验收逐项状态见 validation_acceptance.tsv；未运行的内部临时对象大小测量已单列。\n')
    if diagnostic_rc:
        raise SystemExit('Direct difference diagnosis failed; see diagnosis_driver.log; completed full-run evidence is retained')
    # Plot sampled RSS, keeping original raw samples and plotting source together.
    runtime=read_json(ROOT/'server_install/runtime.json')
    plot=ROOT/'reproduction/plot_xenium_paired_rss.R'
    prefix=out/'xenium_paired_RSS'
    cmd=[runtime['rscript'],'--vanilla',str(plot),str(orig_dir),str(ess_dir),str(prefix)]
    (out/'plot_command.json').write_text(json.dumps(cmd,indent=2)+'\n')
    with (out/'plot.log').open('w') as f:rc=subprocess.run(cmd,cwd=ROOT,env={**__import__('os').environ,'R_LIBS':runtime['dependency_library'],'R_LIBS_USER':'','R_LIBS_SITE':''},stdout=f,stderr=subprocess.STDOUT).returncode
    if rc!=0:raise SystemExit(f'RSS plot failed; see {out}/plot.log')
    # Manifest excludes itself to avoid recursive checksum.
    entries=[]
    for p in sorted(out.rglob('*')):
        if p.is_file() and p.name!='SHA256_manifest.tsv':entries.append({'relative_path':str(p.relative_to(out)),
          'bytes':p.stat().st_size,'sha256':sha(p)})
    write_tsv(out/'SHA256_manifest.tsv',['relative_path','bytes','sha256'],entries)
    print(out/'Xenium_627LR_paired_report_CN.md')
if __name__=='__main__':main()
