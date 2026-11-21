# Figure assets

`FigureS1.{pdf,svg,png}` is the
latest supplemental figure: A–C preserve the reported primary, AutoZyme and
highest-retention results; D uses the completed paired full-627 Xenium RSS curves.
It supersedes the old ESS-only PROFILE panel. The artwork has no visible R1 label.

The final panel A annotation identifies the median of the paired running-time
ratios across 216 records. Panel B identifies the median of the 12
dataset-specific median paired running-time ratios. The final PDF is the
authors' vector master used for both Supplementary Figure S1 and Response
Figure R1.

Run the following from the repository root to redraw all four panels from the
deposited measurements. The Python dependencies are listed in
`supply/environment/plotting_python_versions.txt`.

```sh
python3 supply/figures/rebuild_supplement_figure.py --output s1-redraw
```

The separate, original artwork-finalization step can be replayed on the
historical vector master using Git, Python with pypdf, and Poppler. Run the
following from a checkout that contains commit `67ddb4c`. The PNG and SVG are
exports of the final PDF. The PNG is 600 dpi at A4 size and the SVG uses vector
outlines for the lettering.

```sh
mkdir s1-artwork
git show 67ddb4c1dbe346425c89fc8e0bb535d8b1609be2:supply/figures/FigureS1.pdf > s1-artwork/base.pdf
python3 supply/figures/native_code/finalize_s1_labels.py s1-artwork/base.pdf s1-artwork/FigureS1.pdf
pdftocairo -png -r 600 -singlefile s1-artwork/FigureS1.pdf s1-artwork/FigureS1
pdftocairo -svg s1-artwork/FigureS1.pdf s1-artwork/FigureS1.svg
```

The small paired TSVs are exposed under `../data/` for stable manuscript-facing references.
`response_R1_plot_data/` retains the original response-figure source filenames,
including the independently sampled Original and ESS memory series.

The current Figure 1 is Figure1.pdf. Its final D/E subtitles use the tissue attached to PRJNA289244 (oligodendroglioma). D/E/J retain the same saved outputs, cell-group display mapping and visualization scales. The two-method paired Xenium timing ratios are 13.40 for inference and 8.03 for the whole process. To reconstruct the final main figure, run `python3 supply/figures/main_figure_source/rebuild.py --output /path/to/new-directory` with matplotlib, numpy and PyMuPDF installed. The rendering stages use the archived vector art and original plotting inputs. Original data and base-figure code remain under ../legacy/main_figure_20261101_requires_Xenium_I_update. FigureS1.pdf remains the supplemental figure.
