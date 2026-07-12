# Minimal Reproducible Example

This directory is self-contained and uses a 3.7 MB prepared PBMC3k CellChat object.

- `run_ess.R` runs ESS and checks all outputs against original CellChat.
- `run_three_method_comparison.R` compares original CellChat R, dense Rcpp and ESS.
- `expected_metrics.tsv` defines numerical acceptance thresholds.

From the repository root:

```sh
Rscript examples/minimal_reproducible_example/run_ess.R
Rscript examples/minimal_reproducible_example/run_three_method_comparison.R
```

Both scripts use `nboot = 100` and `seed.use = 1`. Results are written under `results/`.
