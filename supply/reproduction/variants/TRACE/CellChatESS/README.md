# CellChatESS 0.2.2 — release candidate

Efficient Signal Streaming accelerates triMean-based, non-spatial CellChat RNA inference. This source tree is a local candidate prepared on 2026-10-31. It has not been published or deposited as a new release.

## Installation

Install in a private R library with the intended Rcpp and CellChat environment:

```sh
mkdir -p /absolute/path/to/private/library
R CMD INSTALL --preclean --library=/absolute/path/to/private/library /absolute/path/to/CellChatESS
```

```r
.libPaths(c('/absolute/path/to/private/library', .libPaths()))
library(CellChatESS)
# object is a prepared, single-dataset CellChat RNA object.
result <- computeCommunProbESS(object, type = 'triMean', raw.use = TRUE,
    population.size = FALSE, nboot = 100L, seed.use = 1L,
    Kh = 0.5, n = 1, algorithm = 'ess', threads = 1L)
```

The signaling matrix must contain finite, non-negative values and at least one positive value. `algorithm = 'dense'` is the authors' dense C++ comparison path. It is not an independently published competing method. Outputs remain dense: each probability or p-value tensor requires O(n_LR K²) elements.

## Numerical behavior

Single effective ligand/receptor genes now use their group summaries directly, avoiding an unnecessary `exp(log(x))` round trip. Multi-component entities retain the existing geometric mean. Both implementations use the strict `P_boot > P_obs` rule and assign p = 1 when the observed probability is zero. The patch does not change the permutation generator or introduce a tolerance.

Previously recorded targeted tests removed the five known Xenium p-value differences. Those tests used isolated development builds and four implicated LR rows, not a full 627-LR rerun. Probabilities may still differ from CellChat at double-precision rounding scale. Neither exact probability equality for all inputs nor elimination of all possible boundary differences is claimed.

## Version and release status

The revised study uses pinned CellChat 2.2.0.9001 as its common core baseline: original CellChat, the same CellChat with AutoZyme 0.3.1, and this ESS 0.2.2 candidate. The exact common configuration is non-spatial RNA, triMean and population.size=FALSE; this is not a claim of all-feature CellChat compatibility. New primary, AutoZyme, ablation and resource results must be measured from the final source. The original study used CellChat 1.6.1 / ESS 0.1.1; its records are retained as historical supplementary validation. Earlier results labelled 0.2.2 refer to a pre-correction revision build and are tracked separately. Prepared inputs, DB tables and LR sets are shared within each comparison. Source hashes and reproduction instructions are supplied in the surrounding handoff package.

OpenMP is optional. Use `threads = 1L` for the single-thread benchmarks; for parallel checks inspect `object@options$CellChatESS$threads_used` to confirm actual use. A compiler without OpenMP may use one thread even when more are requested.

GitHub: https://github.com/Blake-Deng/CellChatESS

The final release DOI and final validation summaries will be added after server reproduction and author review. Historical source documents are retained in `reference_materials` outside this package.
