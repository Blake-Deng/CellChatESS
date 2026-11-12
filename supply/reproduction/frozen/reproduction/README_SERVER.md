# 已完成记录的重建与独立复现

## 只重建统计/表格/图片（不运行推断）

从本交付根目录运行，Python依赖见environment/plotting_python_versions.txt。主运行证据已解引用复制，读取不依赖服务器符号链接。不要覆盖results/new/core22_statistics。

```sh
python3 reproduction/summarize_and_plot.py --primary results/new/primary216_merged --autozyme results/new/autozyme192 --output /ABS/fresh_statistics
python3 reproduction/assemble_review_tables.py --output .
python3 reproduction/draw_final_revision.py --delivery .
python3 reproduction/write_delivery_documents.py --delivery .
python3 verify_delivery.py
```

已有F_representative_probabilities.tsv时，表格脚本无需R或大RDS。重新导出该概率散点则需从LARGE_FILES.tsv取得两个对应RDS；这些数据来自已保存结果，不需重新计时。表格/图片重建会改变文件哈希，须另存工作副本，保留下载原包用于验证。

## 从源码重新运行（需外部真实输入）

所有入口默认保留已有目录，拒绝覆盖。新机器需要同SHA256的研究输入，不能只凭文件名替代。当前包不含完整研究数据和Xenium大张量。

1. 用environment/server/conda-explicit.txt准备Linux工具链。把所选Conda环境bin放在PATH前。新建私有依赖库，按environment/README_CN.md补齐两项确切版本依赖；保留R包清单和实际安装日志。
2. 创建全新工作根目录，仅复制software、upstream、audit/candidate_source_manifest.tsv、audit/final_candidate.json、reproduction、validation、research_variants（不要复制已完成results/server状态作为新的运行）。不要使用历史已完成状态跳过新任务。
3. `python3 reproduction/install_pinned_sources.py --output /ABS/new_work/server_install --dependency-library /ABS/private_deps --r /ABS/env/bin/R`。该脚本核查三项上游归档及ESS源码包SHA256，输出runtime.json和四个私有库。ESS同一源码在不同R CMD build下的DESCRIPTION元数据可不同，两个tarball均分开保存，生产安装来自software/CellChatESS_0.2.2.tar.gz。
4. 按EXECUTION_PLAN_CN.md运行preflight_versions.R和独立三方法smoke。必须核查实际加载路径、源码body、AutoZyme激活/停用状态和DLL，不能只看DESCRIPTION。原版进程必须保持AutoZyme停用。
5. 在reproduction/config/path_map.example.json的副本填入旧输入前缀→新前缀，将其传给run_unified.py的--path-map。该驱动核验每个输入SHA256。primary216/autozyme192分别使用--limit 1烟测，再分别以--execute和新输出目录完整运行；两批不并行。此次实际最多8对/4输入同时运行，各方法单线程。
6. `python3 reproduction/prepare_variants.py --output server_variants` 从最终生产源码派生A1/A3/A4/TRACE/PROFILE（新工作根目录中此前不能存在server_variants）。按server_campaign.py preparation中的安装、验证、线程探针命令逐一执行；不要沿用旧PASS状态。A2/A5直接使用生产dense/ess路径。
7. 辅助服务器脚本run_server_ablation.py、run_server_traces.py、run_server_xenium.py使用其项目根下server_install/runtime.json；消融清单需做相同输入路径映射。两项Xenium服务器封装保留此次原始输入绝对路径，迁移时只修改该路径并记录diff，固定SHA256。R核心追踪/剖析入口接受输入路径参数，也可直接按实际command.json调用。没有配齐真实输入和受控变体库时不能执行。Xenium profile必须在其他计时队列结束后单独运行。
8. 路径、调度或worker适配需保存纯diff及哈希；数值源码如更改则需新commit/重建/验证。此包未做任何公开发布。

最高保留率输入补测依赖新特征表；run_highest_retention_addendum.py明确选择f_LR最大行，并另存结果。原消融清单rep1不是216记录中最大值，保留为独立历史尝试。重新执行时按此顺序：队列→summarize_and_plot→extract_dataset_characteristics→最高保留率补测；不能拿旧耗时补新结果。

已有直接追踪文件的17位/hex汇合可单独运行 `python3 reproduction/rebuild_trace_exports.py --delivery .`，不运行推断；输出S_G2表，详情见audit/trace_export_note_CN.md。
