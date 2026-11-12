"""Draw paired memory savings as horizontal full-range boxplots on the approved PDF."""
import csv
from pathlib import Path

import numpy as np
import pymupdf as fitz
from matplotlib import font_manager

out = Path(__file__).resolve().parents[1]
baseline = out.parent / 'Figure1_CC22_Xenium_paired_20261111'
input_table = (out.parent / 'Figure1_DEJ_ESS022_bold_20261101' / 'figures' /
               'plot_data' / 'GH_peak_RSS_pairs_192.tsv')
stem = 'Figure1_CC22_paired_H_boxplot_20261111'
with input_table.open() as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
assert len(rows) == 192
assert len({(r['experiment_id'], r['block']) for r in rows}) == len(rows)
for row in rows:
    delta = float(row['autozyme']) - float(row['ess'])
    assert np.isclose(delta, -float(row['ess_minus_autozyme_GiB']))
    assert delta > 0

doc = fitz.open(baseline / 'Figure1_CC22_Xenium_paired_20261111.pdf')
page = doc[0]
state = b'/A1 gs 1 g 0 j 0 w 0 G 1 g\n\n'
start = state + b'244.940642 313.872596 m\n383.075106 313.872596 l'
end = state + b'89.291339 130.110236 m\n304.781102 130.110236 l'
xid, = [x for x in page.get_contents() if start in doc.xref_stream(x)]
stream = doc.xref_stream(xid)
assert stream.count(start) == stream.count(end) == 1
a,b = stream.index(start),stream.index(end)
assert a < b
doc.update_stream(xid,stream[:a]+stream[b:])

layer = fitz.open()
p = layer.new_page(width=page.rect.width, height=page.rect.height)
fontfile = font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans', weight='bold'))
p.insert_font(fontname='PanelBold',fontfile=fontfile)
font = fitz.Font(fontfile=fontfile)
black = (0, 0, 0)
gray = (.3,.3,.3)
ink = (.16,.34,.48)
fill = (.75,.85,.91)
grid = (.88,.90,.91)

def text(x, y, value, size, color=black, align='left'):
    w = font.text_length(value, fontsize=size)
    left = x - {'left':0, 'center':w/2, 'right':w}[align]
    p.insert_text((left,y), value, fontname='PanelBold', fontsize=size, color=color)

text(220.07644,339.9114, 'H', 12.54)
text(244.94064,337.7928, 'Peak RSS saved by ESS', 9.234)
text(244.94064,348.6276, 'ESS below AutoZyme in 192/192 blocks', 6.3, gray)

xmin, xmax = 248., 385.
def xcoord(v):
    return xmin+(xmax-xmin)*v/740

axis_y = 450.
for tick in (0,200,400,600):
    x=xcoord(tick)
    p.draw_line((x,356),(x,axis_y),color=grid,width=.45)
    p.draw_line((x,axis_y),(x,axis_y+2.3),color=black,width=.6)
    text(x,axis_y+10.8,str(tick),6.6,align='center')
p.draw_line((xmin,axis_y),(xmax,axis_y),color=black,width=.7)
text((xmin+xmax)/2,axis_y+24.3,'AutoZyme − ESS peak RSS (MiB)',7.65,align='center')
strata = [('small_1k','1k',72),('medium_10k_or_full','10k or full',72),
          ('large_50k_or_full','50k or full',48)]
for y,(key,label,count) in zip((371.,403.,435.),strata):
    vals=[(float(r['autozyme'])-float(r['ess']))*1024 for r in rows if r['scale']==key]
    assert len(vals)==count
    lo,q1,median,q3,hi=np.quantile(vals,[0,.25,.5,.75,1])
    l,c1,m,c3,h=map(xcoord, (lo,q1,median,q3,hi))
    text(xmin-5,y+2.55,label,7.15,gray,align='right')
    p.draw_line((l,y),(h,y),color=ink,width=.8)
    for v in (l,h):
        p.draw_line((v,y-2.6),(v,y+2.6),color=ink,width=.8)
    p.draw_rect((c1,y-4,c3,y+4),color=ink,fill=fill,width=.65)
    p.draw_line((m,y-4),(m,y+4),color=ink,width=1.15)
    text(m,y-6.9,f'{median:.0f}',7.,ink,align='center')

layer.subset_fonts()
page.show_pdf_page(page.rect,layer,0)
doc.save(out/(stem+'.pdf'),garbage=4,deflate=True)
doc.close()
with fitz.open(out/(stem+'.pdf')) as result:
    pg=result[0]
    pg.get_pixmap(dpi=600,alpha=False).save(out/(stem+'.png'))
    (out/(stem+'.svg')).write_text(pg.get_svg_image(text_as_path=False))
    pg.get_pixmap(dpi=180,alpha=False,clip=fitz.Rect(198,322,411,484)).save(out/'Panel_H_preview.png')
    pg.get_pixmap(dpi=120,alpha=False).save(out/'source/preview.png')

old=(baseline/'Figure1_legend_EN.md').read_text()
first,last=old.index(' (H) '),old.index(' (I) ')
caption=(' (H) Peak process RSS saved by ESS, calculated as the AutoZyme peak RSS minus the ESS peak RSS '
         'within each paired block and expressed in MiB. The boxes span the interquartile range, the vertical '
         'lines mark the medians (labeled to the nearest MiB), and the whiskers extend from the minimum to '
         'the maximum. The 1k, 10k-or-full and 50k-or-full strata contain 72, 72 and 48 paired blocks, '
         'respectively. ESS used less peak memory in all 192 blocks.')
(out/'Figure1_legend_EN.md').write_text(old[:first]+caption+old[last:])
