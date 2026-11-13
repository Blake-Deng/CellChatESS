# Supplementary data and reproducibility deposit

The revision supplement cites the tables and run records in this directory. TSV files retain the recorded precision.
The final CellChat 2.2.0.9001 / ESS 0.2.2 comparisons and historical records are separate.

| Directory | Contents |
|---|---|
| `manuscript/` | Locator for the separate submission documents |
| `data/` | All 216 input characteristics and sampling seeds; 32 AutoZyme inputs; paired results; all 1,008 primary method measurements; dataset/scale summaries and correlations |
| `data/ablation/` | Six-stage medians, all measurements and numerical checks; highest-retention paired repeats |
| `data/validation/` | Structural/input guards, targeted Xenium summaries and all 500 direct decisions in each reference environment |
| `results/` | Compact individual measurement/command/GNU time records, complete numerical comparisons and statistics; full representative method records |
| `xenium/` | Final full-627 paired run, output agreement, memory/time measurements, raw RSS and stage events, frozen scripts and exact identities |
| `environment/` | Recorded Linux Conda explicit lock, R packages/dependencies, sessionInfo, hardware and runtime records |
| `reproduction/` | Frozen scripts, configuration, controlled auxiliary variants and a workspace reconstruction command |
| `figures/` | Current manuscript Figure 1, supplemental Figure S1, native plotting code and source data. |
| `legacy/` | Historical CellChat 1.6.1 / ESS 0.1.1 evidence, distinctly labeled superseded auxiliary records, and the dated main-figure draft requiring panel I synchronization |
| `provenance/` | Source/version identities, copy provenance, SHA256 manifests and explicit exclusions |

## Supplement table map

| Word table | Complete evidence |
|---|---|
| S1 — versions/build scope | `provenance/result_version_build_map.tsv` |
| S2 — dataset characteristics and stratified results | `data/dataset_summary_12.tsv`, `input_characteristics_216.tsv`, `prepared_input_seeds.tsv` |
| S3 — core comparison and uncertainty | `data/primary_pairs_216.tsv`, `autozyme_pairs_192.tsv`, `autozyme_dataset_summary.tsv`, `results/core22_statistics/summary.json` |
| S4 — six-stage ablation | `data/ablation/` |
| S5 — five historical boundary entries | `data/validation/`, `xenium/comparison/five_prior_targets.tsv` |
| S6 — full-627 Xenium pair | `xenium/R2_xenium_memory_time_comparison.tsv`, `RSS_external_stage_summary.tsv`, `audit/theoretical_allocation_payloads.tsv` |
| Figure S1 — supplemental validation | `figures/FigureS1.pdf`, `figures/response_R1_plot_data/` |

The letter's temporary R-table labels do not define this deposit's data filenames.
The Word uses the final S1–S6 numbering. `provenance/copy_provenance.tsv` maps each
copied file to its original handoff path and hash.

The current main figure is `figures/Figure1.pdf`. The accepted submission ledger is `data/Supplementary_Data_S7b.tsv`. The fixed source `data/input_characteristics_216.tsv` contains the same records. For historical preparation of the 216 inputs see `reproduction/input_preparation/`. Complete routine per-run R session information is recoverable from `environment/sessionInfo_by_run/index.tsv` and its three deduplicated contents files.

## Scope of the paired Xenium measurements

Each method's probability tensor and p-value tensor separately contain
3,483,266,523 double elements: 25.95 GiB each, 51.90 GiB combined payload.
The sum refers to the two outputs produced by one method.
Inference clocks cover the function call only. Whole-process clocks additionally
cover startup, input loading, default `aggregateNet`, writing and exit. Peak RSS
is process-wide; in this pair the maximum sampled RSS occurred during aggregation.
Temporary-array sizes estimated from allocation dimensions are theoretical payloads,
not simultaneous measured workspace and must not be summed to reconstruct RSS.
Earlier ESS-only PROFILE values under `legacy/` describe a different scope.

## Verification and reproduction

```sh
python3 supply/verify_supply.py
```

The verification checks every supplied file, all 29 baseline source-manifest entries,
required row counts, mapped evidence paths and exclusions. It does not rerun the
scientific computations. See [`reproduction/README.md`](reproduction/README.md)
before re-execution. The large prepared expression RDS files and full output RDS
objects are not bundled. See [INPUTS.md](../INPUTS.md) for source-data links,
preparation steps and the request route for frozen inputs; manifests retain the expected hashes.
