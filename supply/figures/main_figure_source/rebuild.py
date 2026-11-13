"""Render the accepted main figure from archived vector artwork and measured data."""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from matplotlib import font_manager
import pymupdf as fitz

source = Path(__file__).resolve().parent
supply = source.parents[1]
archived = supply / "legacy/main_figure_20261101_requires_Xenium_I_update"
args = argparse.ArgumentParser(description=__doc__)
args.add_argument("--output", type=Path, required=True, help="A new output directory")
out = args.parse_args().output.resolve()
out.mkdir(parents=True, exist_ok=False)


def stage_copy(src, dest):
    dest.parent.mkdir(exist_ok=True, parents=True)
    shutil.copy2(src, dest)


with tempfile.TemporaryDirectory() as tmp:
    work = Path(tmp)
    base = work / "Figure1_DEJ_ESS022_bold_20261101"
    stage_copy(archived / "Figure1.pdf", base / "output/pdf/Figure1_final_CC22_ESS022_bold_20261101.pdf")
    stage_copy(archived / "plot_data/GH_peak_RSS_pairs_192.tsv", base / "figures/plot_data/GH_peak_RSS_pairs_192.tsv")
    stage_copy(source / "historical_legend.md", base / "figure_legend_EN.md")
    target = Path("CellChatESS_0.2.2_GitHub_with_supply_20261102/supply/xenium/R2_xenium_memory_time_comparison.tsv")
    stage_copy(supply / "xenium/R2_xenium_memory_time_comparison.tsv", work / target)
    stages = [
        ("Figure1_CC22_Xenium_paired_20261111", "update_figure1_i.py"),
        ("Figure1_CC22_paired_H_boxplot_20261111", "update_h.py"),
        ("Figure1_CC22_color_20261111", "update_colors.py"),
    ]
    for folder, script in stages:
        dest = work / folder / "source" / script
        stage_copy(source / script, dest)
        subprocess.run([sys.executable, str(dest)], check=True)
    color = stages[-1][0]
    with fitz.open(work / color / (color + ".pdf")) as doc:
        page = doc[0]
        # The archived artwork supplies the A/F/G labels. Change only their
        # displayed text, alongside the existing D/E tissue-label correction.
        # Panel I's column headings are constructed in update_figure1_i.py.
        spans = [s for b in page.get_text("dict")["blocks"] if "lines" in b
                 for line in b["lines"] for s in line["spans"]]
        fontfile = font_manager.findfont(font_manager.FontProperties(family="DejaVu Sans", weight="bold"))
        font = fitz.Font(fontfile=fontfile)
        replacements = []
        for s in spans:
            x, y = s["origin"]
            size = s["size"]
            if s["text"] == "Glioblastoma · 1k cells · K = 9":
                lines = [(x, y, "Oligodendroglioma · 1k cells · K = 9")]
            elif s["text"] == "Original CellChat vs ESS":
                lines = [(x, 12.8, "CellChat 2.2.0.9001"), (x, 23.25, "vs CellChatESS 0.2.2")]
            elif s["text"] == "Original CellChat":
                lines = [(x, y, "CellChat 2.2.0.9001")]
            elif s["text"] == "Original CellChat probability":
                label = "CellChat 2.2.0.9001 probability"
                lines = [(514.30 - font.text_length(label, fontsize=size) / 2, y, label)]
            else:
                continue
            replacements.append((s, lines))
        if len(replacements) != 6:
            raise ValueError("Expected four archived method labels and two tissue labels")
        for s, _ in replacements:
            page.add_redact_annot(fitz.Rect(s["bbox"]), fill=False)
        page.apply_redactions(images=0, graphics=0)
        for s, lines in replacements:
            for x, y, label in lines:
                page.insert_text((x, y), label, fontsize=s["size"], fontname="ReleaseLabel",
                                 fontfile=fontfile, color=fitz.sRGB_to_pdf(s["color"]))
        doc.save(out / "Figure1.pdf", garbage=4, deflate=True)
