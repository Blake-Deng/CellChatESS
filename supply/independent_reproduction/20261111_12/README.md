# Independent reproduction completed 12 November 2026

[REPORT.md](results/crosscheck/REPORT.md) summarizes the completed replay, and
[final_report.json](results/crosscheck/final_report.json) gives the structured results.
The numerical source, frozen inputs and inference parameters match the manuscript
run. Timings and memory were measured again on the server. Both sets of observed
measurements are reported in the comparison tables.

The directory layout preserves paths relative to the fresh execution root.

- `workspace/results/new/` contains main and AutoZyme campaign completion records and method measurements.
- `results/auxiliary/` contains structural, input-guard, ablation and targeted tracing checks.
- `results/xenium/` contains the full 627-LR comparison and paired resource measurements.
- `results/legacy_cc16_216/` contains the completed legacy-version replay.
- `results/crosscheck/preparation129/` contains input reconstruction and the 216-to-129 ledger.
- `results/figures/` contains figures and statistics generated from this replay.
- `orchestration/compare/` contains the audit and summary scripts used for these records.

All 129 regenerated inputs matched the frozen ordered cells, group assignments,
normalized and signaling expression matrices, databases and LR tables. Their
serialized file hashes differed after resaving. The legacy CellChat 1.6.1 installation
is identified by its recorded environment; its historical upstream Git commit was
not recorded.

The original evidence package has its own manifest. `SHA256SUMS.txt` records the
individual files downloaded from the server summary. The outer supply manifest
also covers the added README, completion note and timing audit.
