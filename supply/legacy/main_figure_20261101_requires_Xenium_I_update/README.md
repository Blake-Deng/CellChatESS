# Dated main Figure 1 draft — panel I update required

This is the 2026-11-01 bold-type main-figure draft, retained with its native plotting
code and source tables. It is **not publication-ready**. Panel I still displays the
earlier instrumented ESS-only PROFILE run (1,105.366 s inference, 1,527.48 s process
time, approximately 111.79 GiB peak RSS). Those measurements do not describe the
final full-627 paired workflow.

The current paired Xenium results are in `../../xenium/`. The current supplemental
figure is `../../figures/FigureS1.pdf`; its panel D uses the completed Nov2 pair.
Do not cite the dated main panel I as the revised resource result or combine its
timings with the new paired run.

A–H/J remain draft artwork assets; D/E/J were generated from version-verified
CellChat 2.2.0.9001 / ESS 0.2.2 output records. The whole main figure must be reviewed
after synchronizing panel I and its caption before publication.

`original_panel_source_map.tsv` preserves the dated source map as provenance, not
the current figure status. `../../figures/asset_status.tsv` records the classification.
The draft native main-figure entry point is `native_code/figure1/redraw_figure1.py`.
It expects a delivery root containing `summary_final.json`, `figures/plot_data/`
and `output/pdf/`; its stored summary remains dated and must not silently be
presented as the new paired evidence. No redraw was performed during this
supplement/repository organization task.
