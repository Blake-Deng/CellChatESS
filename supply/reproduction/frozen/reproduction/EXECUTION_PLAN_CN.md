# 服务器执行计划：四项软件，统一新版核心比较

## 1. 校验、安装、审计

运行 `python3 verify_handoff.py`。软件身份见 `upstream/*_provenance.json`、`audit/final_candidate.json` 和 ESS 的 `inst/source_provenance.json`。三项上游源码及 ESS 均已放入本包；Linux 需要从源码构建，不使用本地 macOS 二进制。

建议在已锁定的 Conda R 环境内，使用新的私有库。`reference_materials/05_code_environment/environment/` 中的锁文件属于历史环境；用于准备依赖时仍需另存本轮实际 Conda explicit、R 包清单、编译器、OpenMP、sessionInfo、CPU/内存与安装日志。不要将历史锁自动标为新运行记录。

```sh
python3 reproduction/install_pinned_sources.py --output /ABS/new/install --dependency-library /ABS/locked/dependency/library --r /ABS/R/bin/R
```

安装器输出四个分开的 R 库和 runtime.json。AutoZyme R 包在其源码归档的 `autozyme_r/` 子目录；CellChat 与 ESS 可用各自源码安装。依赖应按实际锁定环境解决，不能更改上游源码或只改版本号来绕过安装问题。

```sh
R_LIBS=/ABS/locked/dependency/library /ABS/R/bin/Rscript --vanilla reproduction/preflight_versions.R /ABS/new/install/runtime.json /ABS/new/install/extracted/CellChat_2.2.0.9001/CellChat-75253cd0c9e68410e6e721a6d3a0419a1d7e358f /ABS/new/version_audit
python3 validation/run_three_method_smoke.py --runtime /ABS/new/install/runtime.json --input validation/local_smoke/data/pbmc3k_prepared_cellchat.rds --output /ABS/new/server_smoke
```

preflight 检查关键原版函数的实际函数体与上游源码是否一致，再验证 AutoZyme 替换了 computeCommunProb/triMean，停用后恢复。每个 smoke/正式方法运行均为 fresh R 进程。若审计或数值验证失败，保留日志和输入，先诊断再继续；修改 ESS 后提交新源码并重建、重测。

## 2. 论文原有基准思路：在新版基线上重新运行

`primary_216_CC22.tsv` 对应12个数据集、216组既有规模/重复设计（不同prepared输入数量和重复全量输入见清单）。运行 CellChat 2.2.0.9001 原版和最终 ESS 0.2.2，共432个独立方法进程，交替方法顺序。此项重建论文原有速度/内存/数值基准，不是把旧1.6.1结果换标签。

```sh
python3 reproduction/run_unified.py primary216 --runtime /ABS/new/install/runtime.json --output /ABS/new/primary_plan
python3 reproduction/run_unified.py primary216 --runtime /ABS/new/install/runtime.json --output /ABS/new/primary_smoke --limit 1 --jobs 1 --execute
python3 reproduction/run_unified.py primary216 --runtime /ABS/new/install/runtime.json --output /ABS/new/primary216 --jobs 1 --execute
```

先验证全部输入SHA256，路径迁移可加 `--path-map reproduction/config/path_map.example.json` 的本地副本。根据服务器资源安排并发；性能计时避免与其他大计算抢资源，记录实际并发与任务安排。

## 3. 审稿要求的直接三方法比较

`autozyme_32_CC22.tsv` 包含12数据集32份输入。每份执行六种方法顺序，192区组、576个fresh R进程，三种实现为原版 CellChat 2.2.0.9001、该版本激活 AutoZyme 0.3.1、最终 ESS 0.2.2。

```sh
python3 reproduction/run_unified.py autozyme192 --runtime /ABS/new/install/runtime.json --output /ABS/new/three_method_plan
python3 reproduction/run_unified.py autozyme192 --runtime /ABS/new/install/runtime.json --output /ABS/new/three_method_smoke --limit 1 --jobs 1 --execute
python3 reproduction/run_unified.py autozyme192 --runtime /ABS/new/install/runtime.json --output /ABS/new/three_method192 --jobs 1 --execute
```

统一参数：单线程、triMean、raw.use=TRUE、population.size=FALSE、nboot=100、seed.use=1、Kh=.5、n=1。共享同一个prepared对象及DB/LR表，不按方法换数据库。prepared对象可能源自旧预处理；这并不自动成为新版CellChatDB，需如实记录其来源。若新增原生新版预处理或新版数据库输入，应独立生成并记录hash。

计时只分两种：外部包住整个推断函数调用的时间，包含函数内的转换和初始化；GNU time记录整个R进程的墙钟时间与峰值RSS，包含加载、读入、推断、保存与退出。下游过滤/汇总在推断计时之外。默认count比较为prob>0且pval<.05后按群对聚合，不默默增加min.cells过滤。

`compare_results.R` 分别输出概率差异及p/count逐元素比较，p/count不一致会留下表并报告失败。不要为追求完全一致而改成宽容比较。概率的1e-12是小测试接受阈值，不是预设文章中的最大误差；论文报告实际测到的最大差异。

## 4. CellChat 1.6.1：第四项旧版参考

保留 `reference_materials/03_results/benchmark_216` 及对应历史1.6.1/ESS0.1.1环境、源码、日志。新增1.6.1归档锁定到e4f686...，当年运行的准确CellChat commit仍未确认。

允许把历史验证结果作为第四项版本参考单独呈现；不能把历史1.6.1时间与新机器/新运行的另外三项拼成直接性能排名。若需要四项同图比较，另用新库在当前机器、同一输入、线程和计时边界下重跑1.6.1，输出为新运行，并记录其新锁定commit。可用 `validation/smoke_worker.R` 对旧版做小测试；本包本地该项仅为额外兼容检查，不替代历史完整验证。1.6.1下游没有net$LRs字段，ESS写入新版字段属于已记录的接口扩展。

## 5. 结构、输入、支持集与消融

先验证 `validation/run_structure_tests_22.R` 的26项和 `run_input_guards_22.R` 的5类输入，保留合成/真实数据的区别。按 `SUPPORTED_WORKFLOWS.md` 限定支持范围；需要主张旧版兼容时，在单独旧版库重做适用检查。

`extract_dataset_characteristics.R` 计算数据特征、严格正值支持的f_LR与f_union，分母分别为nLR*K²与K²。避免把不同支持集分母混称稀疏度。

`prepare_variants.py` 从最终候选派生A1/A3/A4/PROFILE/TRACE，每个都有纯diff和构建标识；A2为最终algorithm=dense，A5为最终algorithm=ess，A0为原版CC22。A1只取消基因投影，A3取消跨LR缓存和支持剪枝，A4保留缓存但不剪枝。先逐阶段验证输出，再用同一外部函数调用计时。三个代表输入及最高保留率输入见ablation_inputs_CC22.tsv；高保留率A4/A5做三组交替配对，报告实测保留比例和耗时，不预设旧1.78倍。

## 6. 意见6的真实五条目追踪

目标映射/历史五条目的材料在reference_materials及config中。先在历史对应1.6.1环境复核，再在2.2.0.9001环境复核：最终ESS、五个目标各100次，共500次逐次比较。保存两侧置换标签校验、十六进制值、严格判定、计数和p值。`run_target_trace_CC22.R` 默认新版；设置 `CELLCHAT_REFERENCE_VERSION=1.6.1` 并传旧版库可使用旧版。TRACE从最终源码派生，仅插桩，不作性能排名。

本地合成500次追踪只是脚本/构建检查，不能代替实际Xenium。抽取相关LR时保持全体细胞/群并保留全矩阵最大值所在行，保持归一化。五处已知差异消除也不能推断未来所有浮点边界差异都消除。

## 7. Xenium与线程

完整627-LR Xenium输入SHA256为cc6155a5de127ba34b8a26e7f69e7bffc704448f976107253711623d2a88de37。使用最终派生PROFILE构建单线程运行一次，独占足够内存并核存储。每个张量3,483,266,523个double、两个载荷合计约51.90GiB；历史峰值约112GiB仅供资源安排，不是新结果。

`run_xenium_profiled.R` 和独立 `sample_process_memory.py` 记录阶段/RSS；GNU time另记峰值和全进程时间。报告输入对象、两个输出载荷、采样RSS曲线/采样空档、函数时间、全进程时间，并保留稠密O(nLR*K²)限制。当前计划不自动增加原版完整Xenium性能比较；没有严格配对记录时不报告相对原版的加速比。

最终源码在真正支持OpenMP的Linux构建选代表输入做1/4线程数值核查，实际threads_used=4才算4线程。单线程是主要性能比较口径。本地macOS只有串行证据。

## 8. 新结果、图与交付

新图/统计只使用新输出；旧Figure1/R1仅作布局参考。`summarize_and_plot.py`、`export_probability_plot.R` 与 `FIGURE_REBUILD_PLAN_CN.md` 提供起点，最终结合实际数据完成全部正文/补充图表；核对每张图的版本、样本范围、计时范围和数据来源。按数据集抽样计算bootstrap区间，区分pooled与dataset medians。

交付新文件夹/ZIP，含源码/提交、环境锁、全部日志/逐运行表/数值结果、数据与绘图脚本、版本和claim-evidence映射、PASS/FAIL/NOT_RUN/BLOCKED表，以及逐条文稿同步建议。DOI、表号、页码保留给作者，不直接改Word或公开发布。任何新修改都需新commit及对应验证，不能沿用本候选标签掩盖源码变化。
