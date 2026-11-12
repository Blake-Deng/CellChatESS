# Source evidence

CellChat_1.6.1_computeCommunProb.R is a captured function body from the recorded main environment. Line207 applies Pval[Prob==0] <-1 (see zero_probability_source_lines.tsv for machine-extracted exact lines). CellChat_2.2.0.9001_modeling.R is the pinned upstream modeling source. It includes geometric mean, coreceptor positive factors, observed probability construction and the same zero-observed p-value postprocessing. The exact old CellChat installation commit remains unconfirmed; function contents are preserved and hashed.

`original_sampling_seed_rule.py` preserves the seed arithmetic from the original experiment-grid generator; only default output/input paths are made relative. Benchmark seed table is the authoritative dataset-index mapping. This is source evidence, not a command required for packaging.
