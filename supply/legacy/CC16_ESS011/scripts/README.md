# Script use

Run from repository root. Existing benchmark evidence is not overwritten. Optional inference runners require restored `data/` inputs and versioned R libraries; they have not been executed during packaging.

- `recompute_statistics.py`: 10,000 dataset-cluster resamples and feature correlations, reading included rows. Output: results/recomputed. AutoZyme's original 100,000-resample summary is retained with its original separate estimand.
- `run_queues.py`: portable main/compatible queue; plan-only unless `--execute`; main fixed original→ESS with maximum eight pairs; compatible original six order permutations. Each method is a separate process. Compatible full-process measurement wraps Rscript with GNU time; this portable driver does not claim to reproduce historical process startup policy exactly.
- `common_timer_worker.R`: common outer timer, original CellChat and ESS assertions 1.6.1/0.1.1. Input/reference/read/save are outside inference timer. Use `--method --id --input --reference --out --result` pairs.
- `run_compatible_method.R`: five positional arguments METHOD INPUT REFERENCE_OR_DASH OUT_TSV OUT_RDS; asserts compatible versions and one thread.
- `run_ablation_method.R`: METHOD LIBRARY INPUT ORIGINAL_REFERENCE OUTPUT_TSV. A1=`dense_full_gene` on archived ESS0.1.1; A2=`dense_projected` on archived ESS0.1.1; A3=`ess_no_cache` on ESS_no_cache custom library; A4=`ess_no_pruning` on ESS_no_pruning custom library; A5=`ess_pruned` on ESS0.1.1. A0 in the reported table remains saved internal run.time, not a new external wrapper. Each variant library must contain exactly the intended ESS implementation.
- `extract_dataset_characteristics.R`: existing extractor, MANIFEST OUT_TSV SESSION; runs nboot=1 for support counters if invoked. The published characteristic table and observed-only support verification already exist; no invocation is needed for packaging.
- `check_support_sets.R` + `support_observed_only.cpp`: observed-only support masks; excludes bootstrap. Set `CELLCHAT_MAIN_R_LIB`; restore data/inputs and data/timing_outputs. Take a new output directory argument.
- `run_structure_tests_16.R` / `22.R`: existing 26-case drivers; `CELLCHAT_TEST_INPUT` selects the restored 1k input; ESS0.2.2 is required in both cases. Test input and output directory are explicitly separated.
- `run_input_guards_16.R` / `22.R`: INPUT OUTPUT_TSV; ESS0.2.2 and desired CellChat overlay required.
- `run_xenium_profiled.R`: INPUT OUTPUT SUMMARY SESSION, with `CELLCHATESS_STAGE_FILE` pointing to a new stage TSV and `CELLCHAT_ESS_R_LIB` pointing to the instrumented 0.2.2 library. Never load Xenium tensors for plotting; included TSVs suffice.
- `sample_process_memory.py`: `--pid PID --stages NEW_STAGE_TSV --output NEW_RSS_TSV`; samples /proc status/smaps_rollup at one-second intervals; refuses existing output.
- `validate_saved_and_representative.R`: OUTDIR SAVED_TIMER_ROOT REPO_ROOT; P1 reads stored original/ESS arrays, derives default aggregate counts, and P2 reruns only the K9 representative. Restore original figure objects under data/representative. The packaged P2 result is already available.
- `historical_trace/`: source evidence/drivers for five historical boundary comparisons, not a primary revision benchmark. Restore its specific input tensors and instrumented libraries before running; historical complete data are omitted.
- `source_evidence/`: actual CellChat modeling functions, aggregation and original seed-rule source, retained as evidence; not queue entry points.
- `software/`: source-only ESS0.1.1, ESS0.2.2, custom no-pruning/no-cache and profiling snapshots. These are packages/variants, not an assertion that every historical binary had a unique recorded commit.
- `plot_all.sh`: main figure and R1 regeneration, no biological inference.

Package sources and third-party CellChat source evidence remain subject to their original GPL licenses. Source path/library adapters are documented in SOURCE_MANIFEST; they do not change reported numeric results.
