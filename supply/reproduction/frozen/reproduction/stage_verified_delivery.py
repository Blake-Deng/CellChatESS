#!/usr/bin/env python3
"""Copy completed evidence and hash omitted objects; never run inference."""
from pathlib import Path
import argparse,csv,hashlib,json,os,shutil,time,importlib.metadata

ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--delivery',type=Path,required=True);A=p.parse_args()
OUT=A.delivery.resolve();OUT.mkdir(exist_ok=True)
cache={};omitted=[];copied=0
def sha(p):
    st=p.stat();key=(st.st_dev,st.st_ino,st.st_size,st.st_mtime_ns)
    if key not in cache:
        h=hashlib.sha256()
        with p.open('rb') as f:
            for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
        cache[key]=h.hexdigest()
    return cache[key]
def table(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def put(rel,s):
    p=OUT/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(s.strip()+'\n')
def copyfile(src,dst):
    global copied
    dst.parent.mkdir(parents=True,exist_ok=True)
    if dst.exists() and sha(src)!=sha(dst):
        raise RuntimeError('Refusing to replace different delivery evidence: '+str(dst))
    if not dst.exists():shutil.copy2(src,dst)
    copied+=1
def regen(rel):
    if 'primary216' in rel:return 'reproduction/run_unified.py primary216; original experiment_id/method in parent path; same input SHA256'
    if 'autozyme' in rel:return 'reproduction/run_unified.py autozyme192; same input, block and method'
    if 'highest_retention' in rel:return 'reproduction/run_highest_retention_addendum.py; selected_input.tsv'
    if 'ablation' in rel:return 'reproduction/run_server_ablation.py; same stage and prepared input'
    if 'xenium_target_trace' in rel:return 'reproduction/run_server_traces.py; TRACE derivative and five real targets'
    if 'xenium_profile' in rel:return 'reproduction/run_server_xenium.py; PROFILE derivative and full 627-LR input'
    return 'Parent command.json/run.command.json and run.log; retain identical input/source hashes'
def tree(src,rel,skip=(),omit_results=False):
    for root,dirs,files in os.walk(src,followlinks=True):
        dirs[:]=sorted(d for d in dirs if d not in skip and d not in ('__pycache__','.git'))
        for name in sorted(files):
            if name.endswith(('.pyc','.DS_Store')):continue
            q=Path(root)/name;r=Path(rel)/q.relative_to(src)
            if name in skip:continue
            n=q.stat().st_size
            if n>50*1024**2 or (omit_results and q.suffix.lower()=='.rds'):
                omitted.append(dict(delivery_relative_path=str(r),server_path=str(q),bytes=n,sha256=sha(q),reason='bulk output object; all numerical/log records included' if omit_results else 'larger than 50 MiB',regeneration=regen(str(r))))
            else:copyfile(q,OUT/r)

for d in ['upstream','reference_materials','current_author_documents','reproduction','validation','environment','research_variants','logs']:
    tree(ROOT/d,d)
tree(ROOT/'software','software',skip=('CellChatESS_candidate_git',))
# The exact local commit is recoverable from the supplied verified bundle.
tree(ROOT/'server_variants','server_variants',skip=('library',))
tree(ROOT/'server_build_checks','server_build_checks')
tree(ROOT/'server_install','server_install',skip=('library','extracted'))
tree(ROOT/'server_install_attempt01','server_install_attempt01',skip=('library','extracted'))
tree(ROOT/'results/new','results/new',omit_results=True)
for q in (ROOT/'audit').iterdir():
    # New acceptance/claim/version tables are not overwritten by handoff statuses.
    rel=Path('audit/handoff_original')/q.name if q.name in ('acceptance_status.tsv','claim_evidence_map.tsv','result_version_build_map.tsv') else Path('audit')/q.name
    if q.is_dir():tree(q,rel)
    elif q.is_file():copyfile(q,OUT/rel)
for name in ['server_campaign_status.json','server_continuation_status.json']:
    copyfile(ROOT/name,OUT/name)
for name in ['LICENSE','CITATION.cff','.zenodo.json']:
    copyfile(ROOT/'software/CellChatESS'/name,OUT/name)
copyfile(ROOT/'MANIFEST.tsv',OUT/'audit/handoff_original/MANIFEST.tsv')
copyfile(ROOT/'verify_handoff.py',OUT/'audit/handoff_original/verify_handoff.py')
copyfile(ROOT/'server_build_checks/CellChatESS_0.2.2.tar.gz',OUT/'software/server_rebuilt/CellChatESS_0.2.2.tar.gz')
table(OUT/'results/LARGE_FILES.tsv',omitted)
put('results/LARGE_FILES.md',f'''# Bulk objects retained on the server

{len(omitted)} output paths are intentionally omitted, including small per-run RDS outputs as well as large Xenium arrays. See LARGE_FILES.tsv for each path, byte count, SHA256 and regeneration entry. All run logs, measurement tables, comparisons and exact target traces are included. Repeated paths in the merged queue are portable materialized copies of the small evidence; their RDS entries may refer to the same server inode.

Input RDS objects are external research inputs, not software sources. Their expected and actually verified hashes are in audit/server_input_hash_verification.tsv and each queue's input_hash_verification.tsv. Consult reproduction/config and the read-only server provenance for original paths. No source or input file was replaced by packaging.
''')

# Snapshot every installed dependency's file hashes, including non-Conda R packages.
deps=[]
with (ROOT/'environment/server/recursive_dependency_paths.tsv').open() as f:
    for r in csv.DictReader(f,delimiter='\t'):
        folder=Path(r['path'])
        for q in sorted(folder.rglob('*')):
            if q.is_file():deps.append(dict(package=r['package'],version=r['version'],relative_path=str(q.relative_to(folder)),bytes=q.stat().st_size,sha256=sha(q)))
table(OUT/'environment/server/recursive_dependency_file_hashes.tsv',deps)
packages=['numpy','pandas','matplotlib','scipy','tabulate','PyMuPDF','Pillow']
put('environment/plotting_python_versions.txt','\n'.join(f'{n}=={importlib.metadata.version(n)}' for n in packages))
put('reproduction/requirements.txt','numpy>=2.0\npandas>=2.2\nmatplotlib>=3.10\nscipy\ntabulate\nPyMuPDF\nPillow')
put('environment/README_CN.md','''# 环境与构建记录

server/conda-explicit.txt 与 conda-package-list.json 是本轮真实服务器环境，不是历史环境重标。
server/sessionInfo.txt、installed_R_packages.tsv、recursive_dependency_paths.tsv、各方法sessionInfo与DLL哈希构成实际加载证据。
recursive_dependency_file_hashes.tsv 按已安装依赖闭包记录所有文件的SHA256，防止仅凭版本号混同不同构建。

Conda之外，本轮新增私有R依赖 collapse 2.1.8、RcppDist 0.1.1.1；安装日志见 logs/dependencies/install_retry.log。它们不能只靠Conda explicit恢复，需按这些确切版本另装。编译器、Rcpp、BLAS/OpenMP和实际共享库位置均在server及逐运行日志中。共享库哈希仅适用于本次机器/工具链；异机重编译不能预期DLL字节相同。

历史Mac本地验证在validation/local_current与reference_materials中，属于单独环境。实际Linux 4线程证据在results/new/openmp_cc16、openmp_cc22。
''')
put('reproduction/README_SERVER.md','''# 已完成记录的重建与独立复现

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
''')
put('figures/README.md','''# Editable final figure sources

Figure1_final_CC22 and FigureR1_final_CC22 are vector PDF/SVG and 600 dpi PNG. Plot-ready tables are in plot_data/. Code: ../reproduction/draw_final_revision.py. Full captions/alt text: figure_legends.md. Every panel is mapped in panel_source_map.tsv.

The quantitative panels use new CC22/ESS022 measurements. Biological D/E/J remain archived results, with their exact accelerated build unconfirmed; they must not be relabeled as new final-candidate outputs. Their original vectors, group/LR ordering and weights are retained. F is explicitly a newly run different representative input. The Xenium zoom follows new event times (16.5–19 min), not historical timestamps.
''')
put('audit/packaging_note_CN.md',f'''# 本轮交付范围

此包只复制已有实验输出并重建统计/图表，不用旧结果填补新实验。最高保留率addendum是整理时发现清单选错rep后的独立补测，原记录仍保留。

ROOT原始handoff MANIFEST、旧验收表保留在audit/handoff_original；它们的NOT_RUN等状态是交接时状态，不代表服务器最终状态。当前验收表是audit/acceptance_status.tsv。

软件源码Git commit139a61c保持不变；仅复现/整理脚本有新增与适配。research_variants及server_variants不是新的生产版本。生产源码身份见server_source_identity.json、candidate_source_manifest.tsv。全部上游源码归档完整纳入。

本包为作者私有审查交付，原日志保留原始绝对路径以便核对。software/CellChatESS为干净候选，可单独审查；整包还含未发表作者文档和机器信息，不能直接整体公开。未改Word、未填DOI、未进行任何远程发布。

复制小文件计数将在verification报告中核验；大型RDS路径{len(omitted)}条另列，不隐瞒省略。软件和数据的license分别遵循其上游，根GPL-3不覆盖第三方研究输入或作者文稿权利。
''')
print(json.dumps({'copied_small_files':copied,'omitted_output_paths':len(omitted),'dependency_file_hashes':len(deps),'delivery':str(OUT)},indent=2))
