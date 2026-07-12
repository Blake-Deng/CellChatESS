# Methods

## Overview

Expression-Supported Streaming (ESS) is an Rcpp implementation of the triMean-based CellChat communication-probability workflow. It preserves the CellChat object contract and probability equations while changing data preparation, quantile evaluation and supported-pair execution.

## CellChat Quantities

For a gene within one cell group:

```text
triMean(x) = (Q1(x) + 2 * median(x) + Q3(x)) / 4.
```

For ligand summary `L` and receptor summary `R`:

```text
P1 = (L * R)^n / (Kh^n + (L * R)^n).
```

Complex subunits retain geometric-mean aggregation. Co-receptor, agonist and antagonist factors retain the same multiplicative definitions used by CellChat.

## ESS Execution

The signaling matrix remains sparse while ESS computes its maximum and constructs ligand, receptor, complex, cofactor, agonist and antagonist indices. Only referenced rows are selected and converted to the dense matrix required by the C++ kernels. Observed group summaries are calculated only for these projected genes.

Benchmarks compared one full sort with mutable-buffer multi-select. ESS uses one sort below 4096 cells and multi-select at or above 4096 cells. Both implement R type-7 quantile interpolation.

For each interaction, ESS forms ligand-positive sender groups and receptor-positive receiver groups, then evaluates only their Cartesian product. Unsupported entries retain zero probability and p-value one. The probability expression and multiplicative factors are unchanged.

ESS generates bootstrap group labels directly without a separate cell-index permutation matrix. C++ creates bootstrap cell memberships and computes required gene-group triMean summaries on demand. Cached summaries are reused across interactions.

## Complexity Boundary

ESS reduces avoidable internal work, but the returned probability and p-value tensors remain dense with `O(nLR * K^2)` storage. ESS is therefore described as computationally efficient, not as converting the complete CellChat workflow to linear complexity.

## Reproducibility Parameters

Reported validation used `type = "triMean"`, `population.size = FALSE`, `nboot = 100`, `seed.use = 1`, `Kh = 0.5`, `n = 1` and `min.cells = 10`.
