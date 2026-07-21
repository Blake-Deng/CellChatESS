# CellChatESS

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21320633.svg)](https://doi.org/10.5281/zenodo.21320633)
[![Reviewer bundle](https://img.shields.io/badge/reviewer_bundle-v0.2.1-2F5F9E)](https://github.com/Blake-Deng/CellChatESS/releases/download/v0.2.1/CellChatESS_Reviewer_MiniRepo_v0.2.1.zip)

**Efficient Signal Streaming (ESS) for CellChat-compatible communication inference.**

CellChatESS accelerates the triMean communication-probability workflow without changing the CellChat probability model. The repository is intentionally compact: the package source, a self-contained PBMC3k example, and the validation summaries needed to inspect the method are all at the top level.

## What ESS Changes

ESS combines four implementation ideas:

1. Keep `dgCMatrix` input sparse until ligand-receptor, complex and cofactor genes are identified.
2. Compute observed triMean summaries only for genes used by the probability model.
3. Use a benchmark-selected triMean implementation: one sort for groups below 4096 cells and mutable-buffer multi-select for larger groups.
4. Build sender-receiver work from ligand-positive sender groups crossed with receptor-positive receiver groups, then stream bootstrap summaries for supported pairs.

The final CellChat-compatible `K x K x nLR` probability and p-value tensors remain dense. ESS reduces avoidable work and intermediate memory; it does not claim to remove the intrinsic quadratic sender-receiver output term.

## Reviewer Quick Start

For the shortest reproducibility path, download the self-contained
[reviewer mini repository](https://github.com/Blake-Deng/CellChatESS/releases/download/v0.2.1/CellChatESS_Reviewer_MiniRepo_v0.2.1.zip)
and run:

```sh
unzip CellChatESS_Reviewer_MiniRepo_v0.2.1.zip
cd CellChatESS_Reviewer_MiniRepo_v0.2.1
Rscript run_review_demo.R
```

The script checks one- and four-thread execution against the bundled original
CellChat reference and reports a single PASS/FAIL result.

## Installation and compatibility

CellChatESS v0.2.1 was validated with **CellChat 1.6.1**. Install that
CellChat version first from the archived official repository, then install
CellChatESS:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

# Use Bioconductor repositories when resolving CellChat dependencies.
options(repos = BiocManager::repositories())

# The archived upstream commit below contains CellChat 1.6.1.
remotes::install_github(
  "sqjin/CellChat@e4f68625b074247d619c2e488d33970cc531e17c"
)
stopifnot(packageVersion("CellChat") == "1.6.1")

remotes::install_github("Blake-Deng/CellChatESS@v0.2.1")
library(CellChatESS)
```

The validated scope is summarized under [Scope](#scope). Other CellChat
versions or workflow branches should be checked against a reference result
before production use.

Run the repository-integrated v0.2.1 validation after cloning:

```sh
Rscript examples/minimal_reproducible_example/run_review_demo.R
```

This one command checks the serial and conditional-OpenMP paths against each
other and against the bundled original CellChat reference.

## API

```r
cellchat <- computeCommunProbESS(
  cellchat,
  nboot = 100,
  seed.use = 1L,
  threads = 1L,
  algorithm = "ess"
)
cellchat <- computeCommunProbPathwayESS(cellchat)
cellchat <- aggregateNetESS(cellchat)
```

### Thread control

`threads = 1L` uses a true serial fast path and does not enter an OpenMP
parallel region. Values greater than one use OpenMP across bootstrap
permutations within each ligand-receptor interaction:

```r
cellchat <- computeCommunProbESS(
  cellchat,
  nboot = 100,
  seed.use = 1L,
  threads = 4L,
  algorithm = "ess"
)

cellchat@options$CellChatESS[c("threads_requested", "threads_used")]
```

The requested thread count is capped by `nboot`. Builds without OpenMP
automatically use one thread. The dense Rcpp comparator remains serial.

### Benchmark interpretation

CellChatESS supports optional conditional OpenMP beginning with v0.2.0. The
principal runtime and memory results reported in the accompanying manuscript
use **one thread per method** under the same runner and hardware environment.
The four-thread measurements are reported separately as a resource-sensitivity
analysis and are not used to calculate the manuscript's primary speedups.

For the dense Rcpp comparator:

```r
cellchat_dense <- computeCommunProbESS(
  cellchat,
  nboot = 100,
  seed.use = 1L,
  algorithm = "dense"
)
```

## Validation

### 216 matched experiments

Twelve real single-cell datasets were evaluated at six cell-scale settings with three repeats. Each method used `nboot = 100`, `seed.use = 1` and `min.cells = 10`.

| Method | Successful runs | Median core elapsed time | Median peak RSS |
|---|---:|---:|---:|
| Original CellChat R | 216/216 | 426.58 s | 1.861 GiB |
| Dense Rcpp | 216/216 | 35.96 s | 1.159 GiB |
| ESS | 216/216 | 5.18 s | 1.060 GiB |

Core runtime and peak RSS were measured in separate passes. Runtime values
come from the matched core-computation benchmark on common prepared objects;
peak RSS was collected for each method in an independent process using
`/usr/bin/time -v`. Elapsed times observed during the RSS pass are not used in
the runtime summary.

Across all 216 ESS versus original CellChat comparisons:

- maximum absolute probability difference: `1.39e-16`;
- minimum probability Pearson correlation: `1.0`;
- p-values: exact;
- aggregate counts: exact in 216/216 runs;
- pathway and aggregate-weight differences: floating-point precision only.

### Four-thread sensitivity analysis

The conditional-OpenMP implementation was additionally run on all 216 prepared
objects using four native threads and one BLAS thread.

| Threads | Successful runs | Median workflow time | Median peak RSS |
|---:|---:|---:|---:|
| 1 | 216/216 | 5.801 s | 1.060 GiB |
| 4 | 216/216 | 2.660 s | 1.058 GiB |

The paired median 1-to-4-thread speedup was 2.11-fold. Probability, p-value,
pathway, aggregate-weight and count outputs were exactly identical to the
one-thread CellChatESS references for all 216 objects.

### Xenium 50-grid stress test

The largest test contained 1,156,091 cells, 2357 groups and 627 ligand-receptor interactions.

| Method | Core elapsed time | Peak RSS |
|---|---:|---:|
| Original CellChat R | 51,172.58 s | 268.11 GiB |
| ESS | 1,257.72 s | 111.75 GiB |

The complete ESS tensor was compared against the independently validated reference result. Probability, p-value, weight, count and pathway-probability arrays were identical to that reference. Relative to original CellChat, probability differences remained at floating-point precision; five of 3,483,266,523 p-values differed by `0.01`, and none crossed the `0.05` significance threshold.

See [validation/VALIDATION.md](validation/VALIDATION.md) and the TSV files in `validation/`.

## Reviewer mini repository

A self-contained reviewer bundle containing the package source, prepared PBMC3k
object, original CellChat reference, one-command serial/OpenMP check and expected
outputs is attached to the v0.2.1 GitHub Release as
`CellChatESS_Reviewer_MiniRepo_v0.2.1.zip`.

## Repository Map

- `R/`, `src/`: CellChatESS package implementation.
- `examples/minimal_reproducible_example/`: bundled PBMC3k input, scripts and expected metrics.
- `validation/`: concise 216-run and Xenium validation summaries.
- `METHODS.md`: formulas, implementation details and complexity boundaries.
- `CITATION.cff`, `.zenodo.json`: metadata for the independent software release.

## Scope

CellChatESS was validated with CellChat 1.6.1 and currently targets
single-dataset, non-spatial RNA workflows using `type = "triMean"`,
`population.size = FALSE` and `raw.use = TRUE` or a compatible projected
expression matrix.

## Citation

Deng Z. CellChatESS: Efficient Signal Streaming for CellChat inference. Version 0.2.1.
DOI: https://doi.org/10.5281/zenodo.21470360
