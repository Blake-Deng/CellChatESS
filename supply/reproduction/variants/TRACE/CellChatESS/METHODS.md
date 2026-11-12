# Methods

## Overview

Efficient Signal Streaming (ESS) is an Rcpp implementation of the triMean-based CellChat communication-probability workflow. It preserves the CellChat object contract and probability equations while changing data preparation, quantile evaluation and supported-pair execution.

## CellChat Quantities

Let `K` denote the number of cell groups and `n_LR` the number of evaluated
ligand-receptor interactions. The returned probability and p-value arrays have
dimensions `K x K x n_LR`.

For a gene within one cell group:

```text
triMean(x) = (Q1(x) + 2 * median(x) + Q3(x)) / 4.
```

For ligand summary `L` and receptor summary `R`:

```text
P1 = (L * R)^n / (Kh^n + (L * R)^n).
```

A direct ligand/receptor gene uses its group summary directly. Named complexes retain the log-sum geometric mean, even when they contain only one subunit. Measured gene names take precedence over identically named complex rows. Co-receptor, agonist and antagonist factors retain the same multiplicative definitions used by CellChat.

## ESS Execution

Dense numeric matrices and numeric sparse `Matrix` inputs used by CellChat (including `dgCMatrix`) are supported. The signaling matrix remains sparse while ESS computes its maximum and constructs ligand, receptor, complex, cofactor, agonist and antagonist indices. Only referenced rows are selected and converted to the dense matrix required by the C++ kernels. Observed group summaries are calculated only for these projected genes. Other matrix-like object types have not been validated.

Benchmarks compared one full sort with mutable-buffer multi-select. ESS uses one sort below 4096 cells and multi-select at or above 4096 cells. Both implement R type-7 quantile interpolation.

For each interaction, ESS forms ligand-positive sender groups and receptor-positive receiver groups, then evaluates only their Cartesian product. Positive support uses a strict `> 0` criterion. The input expression matrix must be finite and non-negative; negative values, `NA`, `NaN`, and infinite values are rejected before normalization and pruning. Sparse matrix validation examines stored values without converting the matrix to dense form. Under this condition, a ligand, receptor, or required complex subunit with zero group expression has zero aggregate expression and a zero ligand-receptor probability; cofactor, agonist, and antagonist terms multiply that zero and cannot turn it into a positive probability. Unsupported entries retain zero probability and p-value one. Required expression rows must exist; missing complex rows/subunits are rejected. Modifier probability factors use the reference grouping P1 * (agonist_sender * agonist_receiver) * (antagonist_sender * antagonist_receiver).

ESS generates bootstrap group labels directly without a separate cell-index permutation matrix. C++ creates bootstrap cell memberships and computes required gene-group triMean summaries on demand. Cached summaries are reused across interactions.

## Complexity Boundary

ESS reduces avoidable internal work, but the returned probability and p-value tensors remain dense with `O(nLR * K^2)` storage. ESS is therefore described as computationally efficient, not as converting the complete CellChat workflow to linear complexity.

## Reproducibility Parameters

The final-candidate reproduction uses `type = "triMean"`, `raw.use = TRUE`, `population.size = FALSE`, `nboot = 100`, `seed.use = 1`, `Kh = 0.5` and `n = 1`. A downstream minimum-cell filter, where applied, is outside the inference call. The regression count definition is `prob > 0 & pval < 0.05`; no additional minimum-cell filter is silently inserted into that check. Final empirical claims require the server validation records for this exact candidate.
