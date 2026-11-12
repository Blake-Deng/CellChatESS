# Figure assets

`FigureS1.{pdf,svg,png}` is the
latest supplemental figure: A–C preserve the reported primary, AutoZyme and
highest-retention results; D uses the completed paired full-627 Xenium RSS curves.
It supersedes the old ESS-only PROFILE panel. The artwork has no visible R1 label.

The small paired TSVs are exposed under `../data/` for stable manuscript-facing references.
`response_R1_plot_data/` retains the original response-figure source filenames,
including the independently sampled Original and ESS memory series.

The current Figure 1 is Figure1.pdf. Its final D/E subtitles use the tissue attached to PRJNA289244 (oligodendroglioma). D/E/J retain the same saved outputs, cell-group display mapping and visualization scales. The two-method paired Xenium timing ratios are 13.40 for inference and 8.03 for the whole process. To reconstruct the final main figure, run `python3 supply/figures/main_figure_source/rebuild.py --output /path/to/new-directory` with matplotlib, numpy and PyMuPDF installed. The rendering stages use the archived vector art and original plotting inputs. Original data and base-figure code remain under ../legacy/main_figure_20261101_requires_Xenium_I_update. FigureS1.pdf remains the supplemental figure.
