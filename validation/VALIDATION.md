# Validation

## Design

Each paired method used the same prepared CellChat object with `nboot = 100`, `seed.use = 1` and `min.cells = 10`. Peak RSS was measured in an independent process using `/usr/bin/time -v`.

## 216 Matched Experiments

The benchmark comprises 12 real datasets, six target cell scales and three repeats. Original CellChat R, dense Rcpp and ESS completed all 216 runs.

ESS median elapsed time was 5.18 seconds, compared with 38.80 seconds for dense Rcpp and 848.78 seconds for original CellChat R. ESS median peak RSS was 1.060 GiB, compared with 1.159 GiB and 1.861 GiB, respectively.

For ESS versus original CellChat, the maximum probability difference was `1.39e-16`, the minimum Pearson correlation was `1.0`, p-values were exact, aggregate counts were exact in all 216 runs, and pathway/weight differences were limited to floating-point precision.

## Xenium 50-grid Stress Test

The stress test used 1,156,091 cells, 2357 groups, 627 interactions and 100 bootstraps. Original CellChat completed in 51,172.58 seconds with 268.11 GiB peak RSS. ESS completed the core workflow in 1,257.72 seconds with 111.75 GiB peak RSS.

The complete probability, p-value, weight, count and pathway-probability arrays were compared component by component. Relative to original CellChat, probability differences were at floating-point precision. Five of 3,483,266,523 p-values differed by `0.01`; none crossed the `0.05` threshold.

The full large-scale inputs are not redistributed. The bundled PBMC3k example provides a small prepared object and original CellChat reference result.
