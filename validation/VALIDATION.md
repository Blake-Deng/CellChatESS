# Validation

## Design

Each paired method used the same prepared CellChat object with `nboot = 100`, `seed.use = 1` and `min.cells = 10`. Core runtime and peak RSS were collected in separate measurement passes. The runtime summary uses matched core-computation timings; peak RSS was measured for each method in an independent process using `/usr/bin/time -v`. Elapsed times from the RSS pass are not used as runtime results.

## 216 Matched Experiments

The benchmark comprises 12 real datasets, six target cell scales and three repeats. Original CellChat R, dense Rcpp and ESS completed all 216 runs.

ESS median core elapsed time was 5.18 seconds, compared with 35.96 seconds for dense Rcpp and 426.58 seconds for original CellChat R. ESS median peak RSS was 1.060 GiB, compared with 1.159 GiB and 1.861 GiB, respectively.

For ESS versus original CellChat, the maximum probability difference was `1.39e-16`, the minimum Pearson correlation was `1.0`, p-values were exact, aggregate counts were exact in all 216 runs, and pathway/weight differences were limited to floating-point precision.

All principal manuscript comparisons used one thread per method. Conditional
OpenMP was evaluated separately as a resource-sensitivity analysis. Across the
same 216 prepared objects, four-thread CellChatESS completed 216/216 runs with
median workflow time 2.660 seconds and median peak RSS 1.058 GiB, compared with
5.801 seconds and 1.060 GiB at one thread. One- and four-thread probability,
p-value, pathway, aggregate-weight and aggregate-count outputs were exactly
identical in all 216 objects. See `thread_sensitivity_216.tsv`.

## Xenium 50-grid Stress Test

The stress test used 1,156,091 cells, 2357 groups, 627 interactions and 100 bootstraps. Original CellChat completed in 51,172.58 seconds with 268.11 GiB peak RSS. ESS completed the core workflow in 1,257.72 seconds with 111.75 GiB peak RSS.

The complete probability, p-value, weight, count and pathway-probability arrays were compared component by component. Relative to original CellChat, probability differences were at floating-point precision. Five of 3,483,266,523 p-values differed by `0.01`; none crossed the `0.05` threshold.

The full large-scale inputs are not redistributed. The bundled PBMC3k example provides a small prepared object and original CellChat reference result.
