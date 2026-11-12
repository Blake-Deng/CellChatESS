# CellChatESS revision reproduction

Overall status: PASS.

## Computed outcomes

- record_and_array_comparison: exit 0, 0.03 h
- fresh_tables_figures_and_uncertainty: exit 0, 0.01 h
- historical_CC16_ESS011_216: exit 0, 3.61 h
- regenerate_129_prepared_inputs: exit 0, 3.22 h

Core comparison: 3024 same-method arrays checked against historical RDS, 0 changed, 0 historical RDS missing. Within-run numerical acceptance failures: 0.

Legacy 1.6.1/0.1.1: 216/216 pairs completed. Median original/ESS inference ratio: 66.97309481568558 (archived 70.82). The July installed CellChat build is identified by its environment; its upstream Git commit was not recorded.

Prepared inputs: 129/129 regenerated. 129 agreed on ordered cells and groups, normalized and signaling expression, database and LR tables. 0 serialized outputs also matched the frozen SHA256. The 216-to-129 ledger preserves all experimental seeds.

## Re-measured performance

Times, speed ratios and peak RSS are measurements on the current server run. Runtime equality is not required.

- primary / pooled_median: 28.4338 fresh; 27.4771 archived
- primary / median_of_dataset_medians: 23.4899 fresh; 24.542 archived
- autozyme / pooled_median: 1.71751 fresh; 1.75255 archived
- autozyme / median_of_dataset_medians: 1.99637 fresh; 2.0066 archived
- deduplicated_primary / pooled_median: 39.3441 fresh; 38.7593 archived
- deduplicated_primary / median_of_dataset_medians: 35.676 fresh; 36.0298 archived
- ESS faster vs AutoZyme, blocks: 192 fresh; 191 archived
- autozyme median peak RSS, GiB: 1.21944 fresh; 1.22077 archived
- ess median peak RSS, GiB: 0.980024 fresh; 0.981827 archived
- original median peak RSS, GiB: 1.79748 fresh; 1.79671 archived
- Xenium speedup, inference_sec: 10.2891 fresh; 13.4 archived
- Xenium speedup, wall_sec: 6.71703 fresh; 8.03 archived
- Xenium peak RSS reduction percent: 15.4751 fresh; 20.58 archived
- AutoZyme/ESS memory yield count: 192 fresh; 192 archived
- Highest retention A4/A5 median: 1.84022 fresh; 1.83 archived

Full machine-readable statistics, all CIs and individual numerical differences are included in the adjacent tables. Figure1_fresh and FigureS1_fresh were drawn from the new run. The fixed group-display labels and interaction display order in panels D/E/J are reused from the original design.
