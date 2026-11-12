"""Update only panel I of the approved 2026-11-01 vector PDF and its legend."""
import csv
from pathlib import Path

import pymupdf as fitz
from matplotlib import font_manager

out = Path(__file__).resolve().parents[1]
base = out.parent / 'Figure1_DEJ_ESS022_bold_20261101'
evidence = out.parent / 'CellChatESS_0.2.2_GitHub_with_supply_20261102/supply/xenium'
source = base / 'output/pdf/Figure1_final_CC22_ESS022_bold_20261101.pdf'
stem = 'Figure1_CC22_Xenium_paired_20261111'

with (evidence / 'R2_xenium_memory_time_comparison.tsv').open() as f:
    records = {r['method']: r for r in csv.DictReader(f, delimiter='\t')}
orig, ess = records['Original CellChat'], records['CellChatESS']
def values(key):
    return float(orig[key]), float(ess[key])

def comparison(key, digits):
    a, b = values(key)
    return (f'{a:,.{digits}f}', f'{b:,.{digits}f}', f'{a/b:.2f}×')

inference = comparison('inference_function_sec', 3)
wall = comparison('GNU_time_full_process_wall_sec', 2)
rss = values('GNU_time_peak_RSS_GiB')
input_gib = int(orig['input_data_signaling_object_bytes']) / 2**30
prob_gib = int(orig['probability_data_payload_bytes']) / 2**30
pval_gib = int(orig['pvalue_data_payload_bytes']) / 2**30

doc = fitz.open(source)
page = doc[0]
# This dated Matplotlib PDF paints all of I, including its heading, in one
# consecutive stream segment between H's title and the start of J.
stream_id, = page.get_contents()
stream = page.read_contents()
begin = b'1 g 0 j 0.45 w 0.8039215686 0.8235294118 0.8431372549 RG 1 g\n\n445.235615 394.797482 m'
end = b'/A1 gs 1 g 0 j 0 w 0 G 1 g\n\n89.291339 130.110236 m'
assert stream.count(begin) == stream.count(end) == 1
i, j = stream.index(begin), stream.index(end)
assert i < j
doc.update_stream(stream_id, stream[:i] + stream[j:])

# A transparent vector overlay is independently font-subsetted. The other
# panels retain their original graphics streams and embedded font resources.
layer = fitz.open()
p = layer.new_page(width=page.rect.width, height=page.rect.height)
fonts = {}
for weight, tag in [('bold', 'DVB'), ('normal', 'DVR')]:
    fontfile = font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans', weight=weight))
    p.insert_font(fontname=tag, fontfile=fontfile)
    fonts[tag] = fitz.Font(fontfile=fontfile)

black = (0, 0, 0)
gray = (0.32, 0.32, 0.32)
def draw_text(x, y, label, size, weight='DVB', color=black, align='left'):
    width = fonts[weight].text_length(label, fontsize=size)
    if align == 'right':
        x -= width
    elif align == 'center':
        x -= width / 2
    p.insert_text((x, y), label, fontname=weight, fontsize=size, color=color)
    return x, x + width

draw_text(420.3714, 339.9114, 'I', 12.54)
draw_text(445.2356, 337.7928, 'Non-spatial Xenium test', 9.234)
draw_text(445.2356, 348.6276, '1,156,091 cells · K = 2,357', 7.068, color=gray)

# Four columns with units next to each row label. Straight rules and a white
# body keep the two timings visible at the size of the everyday printed figure.
xs = [420.4, 466.6, 512.8, 551.2, 586.4]
ys = [355.6, 376.0, 395.0, 414.0, 433.0]
border = (0.35, 0.36, 0.37)
p.draw_rect((xs[0], ys[0], xs[-1], ys[1]), width=0, color=None, fill=(0.95, 0.95, 0.95))
shape = p.new_shape()
shape.draw_rect((xs[0], ys[0], xs[-1], ys[-1]))
for x in xs[1:-1]:
    shape.draw_line((x, ys[0]), (x, ys[-1]))
for y in ys[1:-1]:
    shape.draw_line((xs[0], y), (xs[-1], y))
shape.finish(color=border, width=0.5)
shape.commit()

def cell(col, row, lines, size, weight='DVB', align='center'):
    top, bottom = ys[row:row+2]
    left, right = xs[col:col+2]
    x = {'left': left+2.8, 'center': (left+right)/2, 'right': right-2.4}[align]
    leading = size+0.6
    baseline = (bottom+top)/2 + size*0.34 - (len(lines)-1)*leading/2
    for line in lines:
        a, b = draw_text(x, baseline, line, size, weight=weight, align=align)
        assert a >= left+0.8 and b <= right-0.8, (row, col, line, a, b)
        baseline += leading

for c, lab in enumerate([['Measure'], ['Original', 'CellChat'], ['ESS'], ['Speedup']]):
    cell(c, 0, lab, 6.6, align='left' if c==0 else 'center')
rows = [(['Inference', 'call (s)'], inference), (['Full-process', 'wall time (s)'], wall),
        (['Peak RSS', '(GiB)'], (f'{rss[0]:.2f}', f'{rss[1]:.2f}', '–'))]
for r, (measure, data) in enumerate(rows, 1):
    cell(0, r, measure, 6.0, align='left')
    cell(1, r, [data[0]], 6.55, 'DVR', 'right')
    cell(2, r, [data[1]], 6.55, 'DVR', 'right')
    cell(3, r, [data[2]], 6.95)
draw_text(xs[0], 446.0, f'Input matrix: {input_gib:.3f} GiB', 7.1)
draw_text(xs[0], 457.0, f'Dense outputs: {prob_gib + pval_gib:.2f} GiB total', 7.1)
draw_text(xs[0], 467.0, f'{prob_gib:.2f} GiB per tensor · 627 LR interactions', 6.65, 'DVR', gray)
layer.subset_fonts()
page.show_pdf_page(page.rect, layer, 0)
meta = doc.metadata
meta.update(title='CellChatESS | Figure 1', subject='Performance and numerical validation, including the paired non-spatial Xenium stress test')
doc.set_metadata(meta)
doc.save(out / (stem+'.pdf'), garbage=4, deflate=True)
doc.close()

with fitz.open(out / (stem+'.pdf')) as final:
    final[0].get_pixmap(dpi=600, alpha=False).save(out / (stem+'.png'))
    (out / (stem+'.svg')).write_text(final[0].get_svg_image(text_as_path=False))
    final[0].get_pixmap(dpi=150, alpha=False).save(out / 'source/preview.png')
    final[0].get_pixmap(matrix=fitz.Matrix(3, 3), clip=(410, 325, 591, 476), alpha=False).save(out / 'source/panel_I_check.png')

old_legend = (base / 'figure_legend_EN.md').read_text()
first = old_legend.index(' (I) ')
last = old_legend.index(' (J) ', first)
legend_i = (
    ' (I) Paired non-spatial Xenium stress test using one single-threaded run per method on the same input '
    '(1,156,091 cells, 2,357 groups, 627 ligand–receptor interactions and 100 bootstrap replicates). '
    'Coordinates were used only to assign cells to a 50 × 50 grid. Spatial inference and distance constraints were disabled. '
    f'Inference time covers the inference function call, including matrix conversion and initialization '
    f'(original CellChat, {inference[0]} s, and ESS, {inference[1]} s, a {inference[2]} speedup). '
    f'Whole-process wall time includes startup, package and input loading, input checks, inference, default network aggregation '
    f'(aggregateNet), saving results and exit ({wall[0]} s and {wall[1]} s, respectively, an {wall[2]} speedup). '
    f'GNU time recorded whole-process wall time and peak resident set size (RSS), which was {rss[0]:.2f} GiB for original '
    f'CellChat and {rss[1]:.2f} GiB for ESS. The input matrix occupied {input_gib:.3f} GiB as an R object. The dense '
    f'probability and p-value tensors each required {prob_gib:.2f} GiB of data storage ({prob_gib+pval_gib:.2f} GiB combined), '
    'excluding object overhead. Storage per tensor scales as O(n_LR K²).'
)
legend = old_legend[:first] + legend_i + old_legend[last:]
legend = legend.replace('# Figure 1 legend - updated biological panels', '# Figure 1 legend', 1)
(out / 'Figure1_legend_EN.md').write_text(legend)
