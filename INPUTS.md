# Inputs for the CellChatESS benchmarks

The package includes a prepared PBMC3k example, the 216-record input ledger,
source accessions, preparation code and the recorded benchmark results. The large
prepared objects used by the main benchmark and Xenium test are obtained
separately.

## Run the included PBMC3k example

Install Rcpp, Matrix and CellChat, then run from the repository root.

```sh
Rscript examples/minimal_reproducible_example/run_review_demo.R
```

The script installs this CellChatESS source into an example-specific library,
computes the original CellChat reference from the supplied input, and checks ESS
with one and four threads requested. It reports the installed package versions
and actual thread counts. OpenMP execution requires a compiler supporting OpenMP.
The prepared input is included at
`examples/minimal_reproducible_example/data/pbmc3k_prepared_cellchat.rds`, with
SHA256 `1305076852f518a6cb57064e8f3678ea3b6c902290206a2463e1d07f2ca817c2`.

## Obtain the frozen benchmark objects

The 216 benchmark records correspond to 129 distinct prepared CellChat objects.
Their identifiers and hashes are in
[`supply/data/Supplementary_Data_S7b.tsv`](supply/data/Supplementary_Data_S7b.tsv)
and
[`supply/reproduction/frozen/reproduction/config/primary_216_CC22.tsv`](supply/reproduction/frozen/reproduction/config/primary_216_CC22.tsv).
The three-method input list is
[`autozyme_32_CC22.tsv`](supply/reproduction/frozen/reproduction/config/autozyme_32_CC22.tsv).

A public bulk-download URL for these large frozen RDS objects is not recorded in
this package. Request the objects through the repository's
[input-request issue page](https://github.com/Blake-Deng/CellChatESS/issues/new).
Include the software version, `dataset_id`, `experiment_id`, target size, repeat
and expected `input_sha256` from the ledger. A request for all 129 distinct
objects can identify the complete `primary_216_CC22.tsv` manifest. Indicate
whether the request also covers the 12 source Seurat objects or the frozen
Xenium object. The private filesystem paths retained in historical manifests
record the original execution locations; use the request page to obtain data.

Verify each supplied prepared object before running the benchmark.

```sh
sha256sum /path/to/prepared_cellchat.rds
```

On macOS, use `shasum -a 256` for this check. The output must match the manifest
entry when replaying the frozen benchmark. Resaving an R object can change its
serialized file hash. Retain the supplied file for the manifest-checked replay.

Restore the recorded script layout and pinned upstream source archives with

```sh
python3 supply/reproduction/make_workspace.py --output /path/to/new-workspace --fetch-upstream
```

This command obtains software sources. Input retrieval uses the request route
above. Configure the Rscript and package-library locations in the reconstructed
workspace's `reproduction/config/runtime.example.json`. Supply a JSON path map
using the original manifest prefix as the key and the local prepared-object
prefix as the value. Then plan and execute a campaign, for example

```sh
python3 /path/to/new-workspace/reproduction/run_unified.py primary216 --runtime /path/to/runtime.json --path-map /path/to/path_map.json --output /path/to/new-results --jobs 1
python3 /path/to/new-workspace/reproduction/run_unified.py primary216 --runtime /path/to/runtime.json --path-map /path/to/path_map.json --output /path/to/new-results --jobs 1 --execute
```

Use `autozyme192` in place of `primary216` for the three-method campaign. Input
hashes, shared parameters, source identities and output comparisons are checked
by the recorded workers.

## Rebuild a prepared input from a source Seurat object

The original preparation script and its 216-record audit are supplied in
[`supply/reproduction/input_preparation/`](supply/reproduction/input_preparation/).
Each audit row identifies the source object, target size, sampling seed and
prepared-object hash. Sampling precedes label selection. Recognized annotation
columns are used when present; otherwise the script normalizes RNA data,
selects variable features, scales, runs PCA, constructs the neighbor graph on
PCs 1–20 and clusters at resolution 0.5. It then applies the human CellChat
database, subsets signaling genes and selects overexpressed genes/interactions.

The recorded preparation environment contains R 4.5.3, Seurat 5.5.1,
SeuratObject 5.4.0, Matrix 1.7-5 and CellChat 1.6.1. Its complete package inventory
is `historical_benchmark_R_packages.txt` in that directory. The original large
Seurat RDS files are external to this package. Their public download URLs and
source-file SHA256 values were not recorded. Obtain the named source objects
through the request route above before using the archival reconstruction code.
The supplied source must include the RNA expression assay, original cell order
and cell metadata used by the preparation script.

For the 10X 1k repeat-1 record, run

```sh
Rscript supply/reproduction/input_preparation/cellchat_benchmark.R --input /path/to/10X.rds --out-dir /path/to/preparation --experiment-id 3CA_data__10X__cells-1000__rep-1__engine-both__ablation-full --dataset-id 3CA_data__10X --engine baseline --n-cells 1000 --seed 112001 --label-col auto --resume true
```

The prepared object is saved under
`/path/to/preparation/checkpoints/3CA_data__10X__cells-1000__rep-1__engine-both__ablation-full/prepared_cellchat.rds`
before that archival script continues with baseline inference and downstream
steps. Use the prepared object and the revised campaign workers for the revised
inference-call timing definition. Use the other audit rows to select each
record's source, target and seed. Request the frozen object when exact file
identity is required and a regenerated serialization differs.

## Source studies and public raw data

The study accessions identify the original data sources. The
[3CA data portal](https://www.weizmann.ac.il/sites/3CA/) supplies curated expression
and annotation files, and the
[3CA methods page](https://www.weizmann.ac.il/sites/3CA/methods) describes its
curation. The portal can change over time. Conversion from currently downloaded
raw matrices to the original 12 Seurat objects requires the original import,
normalization, gene-naming, cell-selection and metadata choices. Those import
scripts and source-object hashes are missing from this deposit. The archived
CellChat preparation code starts from Seurat objects.

| Benchmark dataset | Source object named in preparation audit | Source-study entry |
|---|---|---|
| `3CA_data__10X` | `3CA_data/PRJNA330719/10X/10X.rds` | [Study](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA330719) |
| `3CA_data__EGAS00001004422` | `3CA_data/EGAS00001004422/EGAS00001004422.rds` | [Study](https://ega-archive.org/studies/EGAS00001004422) |
| `3CA_data__Group1` | `3CA_data/PRJNA543322/Group1/Group1.rds` | [Study](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA543322) |
| `3CA_data__Group2` | `3CA_data/PRJNA543322/Group2/Group2.rds` | [Study](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA543322) |
| `3CA_data__PRJNA289244` | `3CA_data/PRJNA289244/PRJNA289244.rds` | [Study](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA289244) |
| `3CA_data__PRJNA395045` | `3CA_data/PRJNA395045/PRJNA395045.rds` | [Study](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA395045) |
| `3CA_data__PRJNA490728` | `3CA_data/PRJNA490728/PRJNA490728.rds` | [Study](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA490728) |
| `3CA_data__PRJNA593551` | `3CA_data/PRJNA593551/PRJNA593551.rds` | [Study](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA593551) |
| `3CA_data__PRJNA638151` | `3CA_data/PRJNA638151/PRJNA638151.rds` | [Study](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA638151) |
| `3CA_data__SmartSeq2` | `3CA_data/PRJNA330719/SmartSeq2/SmartSeq2.rds` | [Study](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA330719) |
| `normal_control__Sample5` | `normal_control/seurat_rds_PRJNA871268/Sample5/Sample5.rds` | [Study](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA871268) |
| `normal_control__Sample8` | `normal_control/seurat_rds_PRJNA871268/Sample8/Sample8.rds` | [Study](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA871268) |

## Xenium input

The source is 10x Genomics'
[Fresh Frozen Human Ovarian Adenocarcinoma with 5K Human Pan Tissue and Pathways Panel](https://www.10xgenomics.com/datasets/xenium-prime-fresh-frozen-human-ovary),
analyzed with Xenium Onboard Analysis 3.0.0. The official page provides output
and supplemental files, including expression and cell-position data, under
CC BY 4.0. Its overview reports 1,157,659 detected cells. The benchmark input
contains 1,156,091 cells after preparation, 648 signaling genes, 2,357 groups
and 627 LR interactions.

The recorded processing normalizes expression to 10,000 counts per cell and
applies log1p, then forms nonempty groups from a 50 × 50 coordinate-quantile grid.
The deposit includes the inference and comparison scripts, parameters and
prepared-input hashes. It does not include the original raw-to-prepared Xenium
construction script or a public download for the exact prepared RDS. Obtain the
frozen Xenium object through the input-request page above for the reported test.
Its expected SHA256 is
`cc6155a5de127ba34b8a26e7f69e7bffc704448f976107253711623d2a88de37`,
as recorded in [`supply/xenium/run_manifest.json`](supply/xenium/run_manifest.json).

The paired inference worker is
`supply/xenium/scripts/run_xenium_paired_worker.R`.
It takes nine positional arguments in this order.

```text
METHOD INPUT INPUT_SHA OUTPUT_DIR CELLCHAT_LIB ESS_LIB DEP_LIB CELLCHAT_MODELING_R ESS_ACCEL_R
```

Use `original` and `ess` as METHOD in separate single-threaded R processes with
the same INPUT and INPUT_SHA. The libraries contain pinned CellChat 2.2.0.9001,
ESS 0.2.2 and their dependencies. The final two arguments identify the pinned
CellChat `modeling.R` and ESS `R/accel_cellchat.R` source files. Set
`EXPECTED_ESS_DLL_SHA256` to the SHA256 of the installed ESS shared library before
the runs. The worker verifies input and source identities, records inference
and aggregation separately, and writes `result.rds`, measurements and stage logs.
Run it under GNU `time -v` to measure whole-process wall time and maximum RSS.
Use `sample_process_resources.py` in the same script directory to collect staged
RSS samples.

After both methods finish, compare their saved arrays with

```sh
Rscript supply/xenium/scripts/compare_xenium_paired_full627.R /path/to/original/result.rds /path/to/ess/result.rds supply/xenium/scripts/xenium_five_targets.tsv /path/to/new-comparison 627
```

Create fresh output directories before invoking the workers and comparator. The
frozen scripts also record their original server paths in commands and manifests;
configure local paths for a new execution. Full input preparation and the
required memory remain substantial, as each probability or p-value array stores
3,483,266,523 double-precision elements.
