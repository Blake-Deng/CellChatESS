"""Adjust marker opacity and the low end of the shared J color scale in the approved PDF."""
from pathlib import Path
import re
import shutil

import numpy as np
import pymupdf
from matplotlib.backends.backend_pdf import pdfRepr
from matplotlib.colors import LinearSegmentedColormap

out = Path(__file__).resolve().parents[1]
base = out.parent / "Figure1_CC22_paired_H_boxplot_20261111"
infile = base / (base.name + ".pdf")
with pymupdf.open(infile) as doc:
    page = doc[0]
    opacity_xref = int(doc.xref_get_key(page.xref, "Resources/ExtGState")[1].split()[0])
    # The baseline plotting code uses A2 for the fills of the individual A/G
    # markers, A6 for the F markers. Check all active uses by their axes clips.
    j_start = b"89.291339 130.110236 m\n304.781102 130.110236 l"
    xid, = [k for k in page.get_contents() if j_start in doc.xref_stream(k)]
    original = doc.xref_stream(xid)
    before_j, j = original.split(j_start, 1)
    before_j += j_start
    allowed = {b"A2": {b"44.6456692913 635.708016599":12, b"44.6456692913 313.8725959423":9},
               b"A6": {b"445.2356147789 468.3938044404":1}}
    for state, clips in allowed.items():
        ops = list(re.finditer(b"/" + state + rb"\s+gs\b", before_j))
        assert len(ops) == sum(clips.values())
        for clip,n in clips.items():
            assert sum(clip in before_j[m.start()-75:m.start()] for m in ops) == n
    for name,old,new in [("A2",0.22,0.36),("A6",0.25,0.39)]:
        assert float(doc.xref_get_key(opacity_xref,name + "/ca")[1]) == old
        doc.xref_set_key(opacity_xref,name + "/ca",str(new))

    old_stops = ["#CDE2FB","#9EC5F4","#6DA7EC","#3987E5","#256ABF","#184F95","#0D366B"]
    new_stops = ["#A9CCF4","#8FBBEE",*old_stops[2:]]
    def pdf_colors(stops):
        cmap = LinearSegmentedColormap.from_list("J",stops, N=256)
        return [tuple(pdfRepr(float(v)) for v in rgb[:3]) for rgb in cmap(np.arange(cmap.N))]
    old_values = pdf_colors(old_stops)
    replacements = dict(zip(old_values,pdf_colors(new_stops)))
    seen = set()
    def update_color(match):
        rgb = match.group(1,2,3)
        assert rgb in replacements
        seen.add(rgb)
        new = replacements[rgb]
        return match.group(0) if rgb == new else b" ".join((*new,match.group(4)))
    # Match both fill and stroke in Matplotlib's two bubble plots and the
    # complete 256-step vector colorbar. PDF operators can span lines.
    color_op = rb"(?<![.\d])([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(rg|RG)\b"
    j = re.sub(color_op,update_color,j)
    assert seen == set(old_values)
    doc.update_stream(xid,before_j+j)
    doc.save(out/(out.name+".pdf"),deflate=True,garbage=4)

with pymupdf.open(out/(out.name+".pdf")) as doc:
    pg=doc[0]
    pg.get_pixmap(dpi=600,alpha=False).save(out/(out.name+".png"))
    pg.get_pixmap(dpi=120,alpha=False).save(out/"source/preview.png")
    (out/(out.name+".svg")).write_text(pg.get_svg_image(text_as_path=False))
shutil.copyfile(base/"Figure1_legend_EN.md",out/"Figure1_legend_EN.md")
