# CellChatESS 0.2.2

Unpublished release candidate, prepared for the CellChat 2.2.0.9001 revision baseline.

- Direct ligand/receptor genes use the group summary directly; named complexes retain the geometric-mean path even when they contain one subunit. Gene names take precedence over identically named complexes, matching the pinned CellChat source.
- Match the probability factor grouping of the reference implementation for agonists and antagonists. Preserve permutation generation and strict bootstrap comparisons without a tolerance.
- Validate expression, cell identities, required LR fields and scalar arguments before native inference. Missing required complex genes are rejected rather than silently omitted; modifier rows use the cofactor table.
- Match the reference pathway output shape when one pathway remains. Document the supported non-spatial RNA/triMean/population.size=FALSE subset.
- Retain optional OpenMP, input gene projection, caching, support-set pruning and streaming. Local smoke results and full server reproduction are reported separately.

# CellChatESS 0.2.1

Documentation and metadata revision of 0.2.0.

# CellChatESS 0.2.0

Optional threads argument and OpenMP bootstrap parallelism.

# CellChatESS 0.1.1

Historical manuscript benchmark snapshot; archived evidence retains its original identity.
