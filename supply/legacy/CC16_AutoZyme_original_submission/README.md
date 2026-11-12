# Historical AutoZyme comparison from the original submission

These three tables are preserved byte-for-byte from the original-submission
reference materials. They describe AutoZyme 0.3.1 and CellChatESS 0.1.1 under
CellChat 1.6.1, with one thread. They are historical records, not the revised
CellChat 2.2.0.9001 / ESS 0.2.2 core comparison.

| File | Historical record |
|---|---|
| `autozyme_216_all.tsv` | 216 AutoZyme runs and saved numerical comparisons |
| `ess_216_rerun_all.tsv` | 216 ESS reruns and saved numerical comparisons |
| `autozyme_vs_ess_peak_216.tsv` | 216 matched AutoZyme/ESS timing and peak-RSS records |

The AutoZyme version is 0.3.1 (commit
`718541d9489596c7c1d75f52e9b3a8b2a429d1f9`). The two numerical-run tables record
CellChat 1.6.1 and CellChatESS 0.1.1 explicitly. AutoZyme was developed for the
newer CellChat environment; these old-environment workflows do not establish its
performance in its supported version. They are not pooled with the revised
192-block comparison in `../../data/autozyme_pairs_192.tsv`. The historical
workflow elapsed times must not be relabeled as the revised inference-call times.
No new computation or timing conversion was performed when preserving these files.

Original local handoff source (relative to the revision workspace):
`CellChatESS_supplementary_agent_input_20261031/02_original_submission_REFERENCE_ONLY/03_supplementary/data/autozyme_headtohead/`.
The original server paths embedded in tables are provenance strings, not portable
execution paths. Copy provenance and exact file SHA256 values are recorded in
`../../provenance/copy_provenance.tsv` and `supply_SHA256_manifest.tsv`.
