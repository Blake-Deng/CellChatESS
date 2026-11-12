# Supported workflows and reference implementation

The revision reference is CellChat 2.2.0.9001, commit
`75253cd0c9e68410e6e721a6d3a0419a1d7e358f`. Compatibility applies to the
subset below; this package is not a complete replacement for all CellChat APIs.

- A preprocessed, single-dataset CellChat object with `options$datatype == "RNA"`,
  `type = "triMean"`, and `population.size = FALSE`.
- Finite, non-negative expression with a positive global maximum; a numeric
  dense matrix or numeric sparse Matrix (including `dgCMatrix`). Stored sparse
  values are validated without densifying the full matrix. Referenced gene
  rows are then projected and converted to the dense native kernel input.
- One non-missing factor identity per expression column, in the same order,
  with no unused levels; unique, non-empty expression gene names.
- `raw.use = TRUE` selects `data.signaling`; `FALSE` selects `data.smooth`.
  Positive integer `nboot`, an integer seed, and positive finite `Kh` and `n`.
- LR rows come from `object@LR$LRsig`. To compare a selected LR subset, prepare
  the identical subset in the input object for all implementations.
- Measured ligand/receptor gene names take precedence over names in the complex
  table. Direct genes use their group summary. Named complexes use geometric
  means, including single-subunit complexes. All required complex genes must
  exist; malformed/incomplete complexes are rejected instead of silently
  deleting missing subunits. Modifiers are resolved through the cofactor table.
- Probability and p-value tensors have dimensions `K x K x nLR`, with source,
  target and interaction names. Bootstrap labels follow the same seeded
  `sample.int()` sequence. Counting remains strictly `P_boot > P_obs`; observed
  zero probabilities have p-value one. No tolerance is applied to comparisons.
- `computeCommunProbPathwayESS()` and the default `aggregateNetESS()` are
  supported. Pathway sums retain p <= threshold; aggregate counts/weights
  retain p < threshold, matching the reference. A single *retained* pathway
  follows the reference dimension-dropping behavior. An input containing only
  one unique pathway can trigger an upstream `aperm()` error; ESS handles that
  input, so this is not claimed as identical upstream error behavior.
- Spatial inference, other mean types, population weighting, merged-object
  inference and additional aggregateNet filtering arguments are outside the
  declared scope. Unsupported inference choices are rejected.

Local smoke checks use real CellChat APIs, a small PBMC input and synthetic
structural/edge fixtures. Exact p/count equality in these checks is an observed
result, not a universal numerical guarantee. Probability values can differ at
floating-point rounding scale. Full final-source benchmarks, known Xenium
boundary tracing and actual OpenMP multi-thread checks belong to server
reproduction. Both output tensors remain dense, requiring O(nLR K^2) storage.
