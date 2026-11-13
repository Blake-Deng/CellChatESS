# Minimal PBMC3k example

This example installs the repository's CellChatESS source into a private library,
computes a fresh reference with the installed CellChat version, and compares
CellChatESS with one thread and with four threads requested. The prepared object
contains fixed expression data, cell labels, database and ligand–receptor tables.

## Run

Install Rcpp, Matrix and CellChat, then run from any working directory.

```sh
Rscript /path/to/CellChatESS/examples/minimal_reproducible_example/run_review_demo.R
```

Use quotation marks around the script path when it contains spaces. A private
`.example_lib/` and a `results/` directory are created beside the script. The
installed CellChat version is reported in the results. The revised core benchmark
uses CellChat 2.2.0.9001. Its pinned source and environment are described in
`../../supply/reproduction/README.md`.

All three inference calls use triMean, raw.use = TRUE, population.size = FALSE,
nboot = 100, seed.use = 1, Kh = 0.5 and n = 1. The script applies the same
minimum-cell filter and the corresponding pathway and network aggregation to
both implementations.

## Checks

The example checks a maximum absolute difference of 1e-12 for probability,
pathway probability and aggregate weight arrays. Empirical p-values and aggregate
counts must agree exactly with the freshly computed original reference. All five
outputs must also be identical between the one- and four-thread ESS requests.

The `threads_used` column records the actual ESS thread count. Compilers without
OpenMP support use the serial fallback.
`results/demo_metrics.tsv`, `results/demo_checks.tsv`, `results/parameters.R` and
`results/sessionInfo.txt` record the run. These example times cover inference and
the downstream steps together and are separate from the manuscript's controlled
inference-call benchmarks.

The prepared input has SHA256
`1305076852f518a6cb57064e8f3678ea3b6c902290206a2463e1d07f2ca817c2`.
Its source is the repository's PBMC3k demonstration prepared in July 2026.
