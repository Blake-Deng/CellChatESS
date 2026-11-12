# 严格配对验证交付索引

本索引仅在队列进入终止状态后生成，不运行推断，不修改原始日志或软件源码。

| 项目 | 状态 | 证据 |
|---|---|---|
| original 全量输出 | COMPLETED；输出哈希核对 True | original/run_status.json |
| ess 全量输出 | COMPLETED；输出哈希核对 True | ess/run_status.json |

完整配对和逐元素比较已完成：是。

完整中文结果：`Xenium_627LR_paired_report_CN.md`。
并列表 R2：`R2_xenium_memory_time_comparison.tsv`；同口径比值：`paired_same_scope_resource_ratios.tsv`。
全部数值、非有限值和原五目标：`comparison/`；差异直接追踪：`difference_diagnostics/`。
RSS 原始采样：各方法 `rss_samples.tsv`；可编辑资源图：`xenium_paired_RSS.*`。

源码与输入身份：`source_and_install_identity.json`、各方法 `run_identity_before.json`；
源码/归档：`sources/`；环境锁：`environment/`；脚本及五目标配置：`scripts/`。
对象大小、纯载荷和 RSS 分开报告；每个稠密输出仍为 O(nLR K²)。
未对干净计时构建的内部临时对象作插桩。源码推导的载荷上界不等于实测同时存活工作区。
终止后的静态文件 SHA256：`terminal_SHA256_manifest.tsv`。动态日志和状态文件另列，不冒称已冻结。
