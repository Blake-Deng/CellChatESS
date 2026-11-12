# CellChatESS 0.2.2

- Support the CellChat 2.2.0.9001 non-spatial RNA workflow using triMean and population.size = FALSE.
- Use cell-group summaries directly for single-gene ligands and receptors, resolving the five diagnosed Xenium p-value differences. Preserve CellChat's permutation sequence and strict bootstrap comparison rule.
- Align complex, cofactor, agonist, antagonist and pathway calculations with the pinned CellChat reference.
- Validate finite nonnegative expression input, cell identities, LR fields and scalar arguments. Check sparse expression matrices without densification.
- Retain optional OpenMP parallelism and the serial inference path.
- Provide manuscript results, pinned environments, source identities and a complete independent replay covering 216 main pairs, 192 three-method blocks, 627 Xenium LR interactions, 216 legacy pairs and reconstruction of 129 distinct inputs.

# CellChatESS 0.2.1

Documentation and metadata revision of 0.2.0.

# CellChatESS 0.2.0

Optional threads argument and OpenMP bootstrap parallelism.

# CellChatESS 0.1.1

Original manuscript benchmark version.
