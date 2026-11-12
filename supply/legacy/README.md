# Historical and superseded evidence

`CC16_ESS011/` preserves historical CellChat 1.6.1 / ESS 0.1.1 verification.
Its version-specific results do not define the revised CellChat 2.2.0.9001 core comparison.
The exact original CellChat 1.6.1 source commit was not recoverable from historical records;
the separately pinned 1.6.1 source in `provenance/upstream/` is a new frozen reference.

`CC16_AutoZyme_original_submission/` preserves the original-submission AutoZyme
0.3.1 / ESS 0.1.1 comparison under CellChat 1.6.1. Its three tables are historical
reference records, separate from the compatible-environment revised comparison;
see its README for timing and compatibility limits.

`superseded_ESS_only_PROFILE/` is an earlier instrumented ESS-only resource run.
Its 1,105.366 s inference, 1,527.48 s process time and 111.79 GiB peak must not be
combined with the final full-627 paired workflow in `../xenium/`.
`superseded_version_map/` and `earlier_retention_rep1/` are retained solely for provenance.
The current version map is `../provenance/result_version_build_map.tsv`.
