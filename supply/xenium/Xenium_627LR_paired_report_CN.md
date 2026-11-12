# Xenium 全量 627 LR 配对验证报告（服务器复现）

本报告只针对下述冻结 Xenium 输入、源码构建、服务器和参数，不外推到任意输入。原版 CellChat 与最终 ESS 在独立 R 进程中依次运行，未并发。

## 输入、软件与参数

- 输入：`/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/prepared_xenium_ovary_full_grid50.rds`；SHA256 `cc6155a5de127ba34b8a26e7f69e7bffc704448f976107253711623d2a88de37`。两次运行在读入前和运行中核对相同输入。
- 对象：1,156,091 细胞、648 信号基因、K=2,357、627 LR；两次运行的 cell/group/DB/LR 序列化哈希完全一致。
- Original：CellChat 2.2.0.9001；上游 commit `75253cd0c9e68410e6e721a6d3a0419a1d7e358f`；归档 SHA256 `52a90ce9f121acd58c9e5b2542956de0f3d81ac9b10906454847946f289e0739`。实际 `computeCommunProb` 函数体与冻结源码一致，AutoZyme 未加载。
- ESS：CellChatESS 0.2.2；commit `139a61c33d4c09568f48f4afd0a16dc57bb22f4f`；tree `62fb33d233e8b2296eca11b7eefade7745429e10`；安装 DLL SHA256 `c46577016d28829e5dcd1c694b3f185a762a6f925bf05cfe0721981030085752`。实际 `computeCommunProbESS` 函数体与冻结源码一致。
- 参数：RNA、`type="triMean"`、`raw.use=TRUE`、`population.size=FALSE`、`nboot=100`、`seed.use=1`、`Kh=0.5`、`n=1`；Original 与 ESS 均单线程，ESS `algorithm="ess"`。

## 时间与内存

- 推断函数时间从各函数调用前计时，函数返回后停止；不含 R 启动、包加载、输入读取、`aggregateNet` 或序列化。
- GNU time 全进程墙钟从 R 进程启动到退出，包含推断函数时间、同一 `aggregateNet` 默认计数步骤和相同 gzip `saveRDS`；全进程时间包含推断时间，不能再相加。
- RSS 同时报告 GNU time 的进程峰值和独立 1 秒 `/proc` 原始采样峰值。采样间隔与大于 2/5 秒的缺口见 `rss_sampling_gap_summary.tsv`，没有插值。
- `object.size`、稠密数组纯 payload 和进程 RSS 是不同量。每个概率/p 值输出仍是稠密 `n_LR × K × K` 数组，空间复杂度为 O(n_LR K²)；ESS 支持集剪枝不消除这两个稠密输出。

## 全量逐元素核查

- 概率与 p 值维度、名称和顺序一致：是 / 是。
- 概率元素数 3,483,266,523；数值不同元素 2,318,979；最大绝对差 `2.2204460492503131e-16`。最大差坐标：LR `THBS1_CD47`，发送 `Group_Grid_50_24`，接收 `Group_Grid_22_37`；Original `0.50356041579253241`（`0x1.01d2abbac736ap-1`），ESS `0.50356041579253263`（`0x1.01d2abbac736cp-1`）。
- p 值元素数 3,483,266,523；不一致元素 0。所有差异位置及十进制/十六进制值见 `comparison/pvalue_differences.tsv`。未使用容差。
- 默认 `aggregateNet(thresh=0.05)` count：Original 和 ESS 的独立分块重算均与各自保存值吻合：是 / 是；两方法 count 矩阵完全一致：是。count 差异坐标见 `comparison/aggregate_count_differences.tsv`。
- 概率及 p 值的 NA、NaN、正负 Inf 计数见 `comparison/nonfinite_summary.tsv`；未有限值差异坐标见对应文件。

## 原五个目标条目

按 LR、发送群、接收群三者共同定位，完整输出数值和 double 十六进制值见 `comparison/five_prior_targets.tsv`。
- CSF1_CSF1R / Group_Grid_3_32 → Group_Grid_12_50：P Original=0.0066836308070305859，ESS=0.0066836308070305859；p Original=0.55000000000000004，ESS=0.55000000000000004；count=4/4。
- LAMA3_CD44 / Group_Grid_17_26 → Group_Grid_2_45：P Original=0.051898326898552324，ESS=0.051898326898552324；p Original=0.59999999999999998，ESS=0.59999999999999998；count=1/1。
- LAMA3_CD44 / Group_Grid_25_16 → Group_Grid_2_45：P Original=0.05324220677168226，ESS=0.05324220677168226；p Original=0.56000000000000005，ESS=0.56000000000000005；count=2/2。
- LAMC1_CD44 / Group_Grid_12_50 → Group_Grid_10_19：P Original=0.0064904606955402703，ESS=0.0064904606955402703；p Original=0.81999999999999995，ESS=0.81999999999999995；count=11/11。
- LAMC2_CD44 / Group_Grid_12_50 → Group_Grid_10_19：P Original=0.0064904606955402703，ESS=0.0064904606955402703；p Original=0.56999999999999995，ESS=0.56999999999999995；count=11/11。

全量 627 LR 的经验 p 值无差异；因此本次全量结果未触发追加 bootstrap 逐次追踪。原五目标此前的目标级追踪记录仍是单独范围证据。

## 运行文件

- 每种方法的原始 R 日志、GNU time、RSS 采样、函数源码快照、sessionInfo、输入对象 group sizes/LR 表和 `run_metadata.json` 位于 `original/` 与 `ess/`。
- R2 并列数据位于 `R2_xenium_memory_time_comparison.tsv`；RSS 绘图源数据与绘图脚本见 `*_rss_plot_data.tsv` 和 `reproduction/plot_xenium_paired_rss.R`。
- 本次比较脚本：`reproduction/compare_xenium_paired_full627.R`；worker：`reproduction/run_xenium_paired_worker.R`；队列入口：`reproduction/run_xenium_paired_full627.py`。
- 旧的 ESS PROFILE 插桩运行没有被用于本次 ESS 数值结果或配对时间/RSS。

## 原始输出和复现边界

完整 RDS 输出在本目录的 `original/result.rds` 与 `ess/result.rds`。所有结论仅适用于该冻结输入、上述软件源码、硬件和参数；概率逐元素相同以外的结论需按输出报告。

## 本次实测时间与 RSS

| 方法 | 函数调用秒 | 全进程秒 | GNU time 峰值 GiB | 采样峰值 GiB |
|---|---:|---:|---:|---:|
| Original CellChat | 17186.512 | 18092.000 | 250.07668 | 250.12001 |
| CellChatESS | 1282.588 | 2251.840 | 198.61917 | 198.62785 |

同口径推断函数耗时比 13.399870；同口径全进程耗时比 8.034319；GNU time 峰值 RSS 降幅 20.576692%。只针对本次各一次顺序运行，不提供重复实验置信区间。

主性能构建未对函数内部临时对象作插桩，外部阶段 RSS 见 RSS_external_stage_summary.tsv。辅助 PROFILE 构建的阶段记录单独保留，不能替代本次主构建的时间/RSS，也不能把峰值 RSS 减去输出载荷解释为测得的临时工作区。独立 RSS 采样无法观察进程启动前或退出后的内存，采样缺口不会以插值补齐。

## 聚合 count 的两个明确口径

默认 aggregateNet 显著性计数与原稿 filterCommunication(min.cells=10) 后的计数分别核查。官方 CellChat modeling.R:916-920 排除的是细胞数小于或等于 10 的群；min.samples=NULL 的有效值为 1，因此不触发多样本过滤。
本次原稿口径排除 24 个群对应的发送/接收页；保留完整 K×K 矩阵并将这些行列置零。Original 显著 LR 边总数为 69890010，ESS 为 69890010；群对计数差异数 0，逐元素一致：是。
细胞数过滤的独立核查在主性能计时之外；GNU time 全进程范围仍是共同 worker 中的输入、推断、默认未加细胞数过滤的聚合、保存与退出，不能将该全进程时间称为包含细胞数过滤的完整论文工作流。规则、排除群列表和两份计数矩阵见 comparison/aggregate_count_policy_summary.tsv 及相邻文件。

## 差异追踪与验收

差异诊断状态：NOT_REQUIRED。若存在差异，逐 LR 的直接 bootstrap 追踪与全量输出交叉核对保存在 difference_diagnostics；目标级诊断不替代全量比较，也不用于计时。具体浮点来源需要根据直接追踪审阅，不预设为 BLAS、编译器或求和顺序。
验收逐项状态见 validation_acceptance.tsv；未运行的内部临时对象大小测量已单列。
