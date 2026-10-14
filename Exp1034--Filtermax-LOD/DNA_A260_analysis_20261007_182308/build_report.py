"""Create a printable lab-notebook PDF and browsable HTML from the audited analysis."""
from html import escape
import json
from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak
from PIL import Image as PILImage

HERE = Path(__file__).resolve().parent
FONT = Path('/usr/share/fonts/truetype/dejavu')
pdfmetrics.registerFont(TTFont('DejaVu', str(FONT / 'DejaVuSans.ttf')))
pdfmetrics.registerFont(TTFont('DejaVu-Bold', str(FONT / 'DejaVuSans-Bold.ttf')))
pdfmetrics.registerFontFamily('DejaVu', normal='DejaVu', bold='DejaVu-Bold')
INK = colors.HexColor('#183341')
TEAL = colors.HexColor('#007E75')
STYLES = {
  'title': ParagraphStyle('title', fontName='DejaVu-Bold', fontSize=24, leading=29, textColor=INK, spaceAfter=13),
  'h2': ParagraphStyle('h2', fontName='DejaVu-Bold', fontSize=15, leading=19, textColor=TEAL, spaceBefore=8, spaceAfter=8),
  'body': ParagraphStyle('body', fontName='DejaVu', fontSize=10.5, leading=15, textColor=INK, spaceAfter=9),
  'small': ParagraphStyle('small', fontName='DejaVu', fontSize=9, leading=12.5, textColor=INK, spaceAfter=7),
  'table': ParagraphStyle('table', fontName='DejaVu', fontSize=9, leading=12, textColor=INK),
}


def build(root=HERE):
  meta = json.loads((root / 'analysis_metadata.json').read_text())
  figures = json.loads((root / 'figure_index.json').read_text())
  g = pd.read_csv(root / 'tables/group_summary.csv')
  b = pd.read_csv(root / 'tables/blank_statistics.csv')
  fits = pd.read_csv(root / 'tables/linear_fits.csv')
  dep = pd.read_csv(root / 'tables/group_departures.csv')
  limits = pd.read_csv(root / 'tables/exploratory_limits.csv')
  story, html = [], []
  width, height = landscape(A4)
  usable = width - 80

  def paragraph(text, style='body'):
    story.append(Paragraph(text, STYLES[style]))
    tag = 'h1' if style == 'title' else 'h2' if style == 'h2' else 'p'
    html.append(f'<{tag}>{text}</{tag}>')

  def table(headers, rows, widths=None):
    values = [[Paragraph(escape(str(x)), STYLES['table']) for x in row] for row in [headers]+rows]
    t = Table(values, colWidths=widths or [usable/len(headers)]*len(headers), repeatRows=1, hAlign='LEFT')
    t.setStyle(TableStyle([
      ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#DCEEEB')),
      ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.white, colors.HexColor('#F3F6F7')]),
      ('VALIGN', (0,0), (-1,-1), 'TOP'), ('TOPPADDING', (0,0), (-1,-1), 7),
      ('BOTTOMPADDING', (0,0), (-1,-1), 7), ('LINEBELOW', (0,0), (-1,0), .6, TEAL),
    ]))
    story.extend([t, Spacer(1, 10)])
    html.append('<table><thead><tr>'+''.join(f'<th>{escape(str(x))}</th>' for x in headers)+'</tr></thead><tbody>'+
                ''.join('<tr>'+''.join(f'<td>{escape(str(x))}</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table>')

  def page():
    story.append(PageBreak())
    html.append('<hr>')

  paragraph('Exp1034 | DNA A260 volume & sensitivity', 'title')
  paragraph('Lab notebook analysis • 7 October 2026 • FilterMax F5, serial 1191', 'small')
  paragraph('<b>Practical conclusion:</b> Start with <b>100 µL</b> for low-concentration work. The lowest tested concentration, '
            '<b>2.1875 ng/µL</b>, produced a signal above the matched blanks and an 8.5% triplicate CV at 100 µL. '
            'Doubling volume to 200 µL gave 8.0% CV and almost twice the signal, but also higher blank variation.')
  paragraph('Concentration correction used throughout', 'h2')
  paragraph('The user clarified that A30 contained <b>70 ng/µL</b>. The exported metadata still records 16 ng/µL. '
            'This analysis preserves those original fields and uses the corrected nominal twofold series below. '
            'Concentration slopes decrease by 4.375-fold and concentration-limit estimates increase by 4.375-fold; '
            'blank-corrected readings, CVs, R² and volume comparisons are unchanged.')
  table(['Source', 'A30', 'A29', 'A28', 'A27', 'A26', 'A25', 'A31'],
        [['Nominal ng/µL', '70', '35', '17.5', '8.75', '4.375', '2.1875', '0 (blank)']])
  paragraph('Lowest tested level: what happened at each volume?', 'h2')
  rows=[]
  interpretation={50:'Inconsistent; D9 elevated',75:'Detectable, less precise',100:'Best supported starting choice',200:'More signal; little precision gain'}
  for r in g[g.concentration_ng_uL == 2.1875].itertuples():
    rows.append([f'{r.volume_uL:g}', f'{r.mean_corrected:.5f}', f'{r.sample_sd:.5f}', f'{r.cv_percent:.1f}%', interpretation[r.volume_uL]])
  table(['Volume (µL)', 'Corrected OD', 'Sample SD (OD)', 'CV, n=3', 'Interpretation'], rows,
        [90,110,120,90,usable-410])
  paragraph('All 96 wells are present and finite: 72 DNA wells and 24 A31 blanks. No wells are excluded. '
            'Electronic sample IDs, sources, volumes and replicates match the saved layout. This does not independently verify physical contents.', 'small')
  page()

  paragraph('How low can the volume and concentration go?', 'title')
  paragraph('<b>50 µL:</b> the minimum tested volume. Its blanks are quiet, and concentrations 8.75-70 ng/µL have CVs below 6.4%. '
            'At 2.1875 ng/µL, however, D9 raises the triplicate CV to 100.5%. Do not generalize high-concentration performance to the lowest level. '
            '<b>75 µL:</b> a volume-saving candidate; CV is below 4% at 8.75-70 ng/µL, but 18.9% at 2.1875 and 12.9% at 4.375 ng/µL.')
  paragraph('Corning lists a recommended working range of 75-200 µL for plate 3635 [1]. Thus 50 µL is an exploratory use below that range; '
            'it is retained here, with its observed limitations. Neither this specification nor this single plate proves a universal minimum volume.')
  paragraph('Detection is demonstrated only at the tested levels', 'h2')
  paragraph('At 2.1875 ng/µL, all three wells at every volume exceed the matching blank mean plus 3.3 blank SD. '
            'This is evidence of separation from these blanks, not a validated detection probability. No DNA sample below 2.1875 ng/µL was tested. '
            'The noise-based estimates below are targets for follow-up, not demonstrated limits.')
  rows=[]
  for v in [50,75,100,200]:
    t=limits[(limits.volume_uL==v)&(limits.fit_kind=='concentration_full')].set_index('metric')
    fmt=lambda metric:f'{t.loc[metric,"estimate_ng_uL"]:.2f} ({t.loc[metric,"ci95_low"]:.2f}-{t.loc[metric,"ci95_high"]:.2f})'
    rows.append([v,fmt('DL'),fmt('QL'),'8.75' if v in [50,75] else '4.375'])
  table(['Volume (µL)', 'Noise DL (95% interval)', 'Noise QL (95% interval)', 'Lowest tested passing strict screen'],rows)
  paragraph('All concentration values above are ng/µL. DL = 3.3 × blank SD / slope; QL = 10 × blank SD / slope [2]. '
            'Intervals resample wells on this plate only. The illustrative, post hoc <b>strict screen</b> requires all three wells above 10 blank SD '
            'and a corrected CV ≤10%; it is not a validated LOQ or a predefined acceptance criterion.', 'small')
  paragraph('For perspective, 2.1875 ng/µL uses 109.375, 164.063, 218.75 and 437.5 ng at 50, 75, 100 and 200 µL. '
            'These masses assume the stated concentrations and planned volumes. No path-length correction was applied; absorbance cannot be converted '
            'to an independent DNA concentration or purity result from this read alone [3].', 'small')
  page()

  for figure in figures:
    image_path = root / figure['png']
    iw, ih = PILImage.open(image_path).size
    scale=min(usable/iw, 417/ih)
    story.append(Image(str(image_path), width=iw*scale, height=ih*scale, hAlign='CENTER'))
    story.append(Spacer(1, 7))
    story.append(Paragraph(figure['caption'], STYLES['body']))
    html.append(f'<section><img src="{figure["png"]}" alt="{escape(figure["title"])}"><p>{escape(figure["caption"])}</p>'
                f'<p class="download"><a href="{figure["png"]}">Download PNG</a> · <a href="{figure["svg"]}">Download SVG</a></p></section>')
    page()

  paragraph('Regression results | free intercept throughout', 'title')
  paragraph('Concentration fits use nominal ng/µL on the x-axis; slopes are OD/(ng/µL). Full range: 2.1875-70, n=18 wells per fit. '
            'Lower four levels: 2.1875-17.5, n=12. Blank wells are not fitted as calibration points. RMSE = √(SSE/n).', 'small')
  rows=[]
  for r in fits[fits.fit_kind!='volume'].itertuples():
    rows.append([f'{r.volume_uL:g}', 'Full' if r.fit_kind=='concentration_full' else 'Lower four',
                 f'{r.slope:.6f}',f'{r.slope_ci95_low:.6f}-{r.slope_ci95_high:.6f}',
                 f'{r.intercept:.5f}',f'{r.r_squared:.5f}',f'{r.residual_rmse:.5f}'])
  table(['µL','Range','Slope','Slope 95% interval','Intercept','R²','RMSE (OD)'],rows,[50,90,85,195,85,75,usable-580])
  paragraph('All slope intervals are above zero. The full curves are strongly increasing and approximately straight overall. '
            'Nevertheless, systematic group departures are visible at 75, 100 and 200 µL. At 50 µL the low-range fit deteriorates sharply '
            '(R²=0.726), largely because of D9. A high full-range R² is not enough to justify low-end quantitation.')
  paragraph('Linearity conclusions use the group-residual intervals on the next page, not an arbitrary R² cutoff. '
            'No simultaneous or multiple-comparison adjustment is applied; these departures are exploratory.', 'small')
  page()

  paragraph('Volume fits and where straight lines depart', 'title')
  paragraph('Volume fits use 50, 75, 100 and 200 µL at each concentration, n=12 wells. Slope units are OD/µL.', 'small')
  rows=[]
  for r in fits[fits.fit_kind=='volume'].itertuples():
    q=dep[(dep.fit_kind=='volume')&(dep.concentration_ng_uL==r.concentration_ng_uL)&dep.excludes_zero]
    flag=', '.join(f'{v:g}' for v in q.volume_uL)+' µL' if len(q) else 'No clear departure detected'
    rows.append([f'{r.concentration_ng_uL:g}',f'{r.slope:.7f}',f'{r.intercept:.5f}',f'{r.r_squared:.4f}',f'{r.residual_rmse:.5f}',flag])
  table(['ng/µL','Slope','Intercept','R²','RMSE (OD)','Group intervals excluding zero'],rows,[75,100,95,65,95,usable-430])
  paragraph('Concentration fits: groups driving departures', 'h2')
  statements=[]
  for v in [50,75,100,200]:
    q=dep[(dep.fit_kind=='concentration_full')&(dep.volume_uL==v)&dep.excludes_zero]
    items=', '.join(f'{r.concentration_ng_uL:g} ({"above" if r.departure>0 else "below"})' for r in q.itertuples())
    statements.append(f'<b>{v} µL:</b> {items} ng/µL')
  paragraph('; '.join(statements)+'. These are group means whose conditional 95% departure intervals exclude zero. '
            'The direction is relative to the fitted line, not to a known true concentration.', 'small')
  paragraph('The 2.1875 ng/µL volume fit has R²=0.321 and wide uncertainty because the 50 µL group is inconsistent. '
            'Its lack of a clear group departure is not affirmative evidence of linearity. At higher concentrations, the small departures may reflect '
            'optics, dilution accuracy, dispensing differences or individual wells; this plate does not distinguish the causes.', 'small')
  page()

  paragraph('Methods, audit trail and next experiment', 'title')
  paragraph('<b>Acquisition:</b> run dna-a260-f5-20261007T182308Z-4a000d5d; 260 nm, 12 nm bandwidth, excitation slide 2/slot 1; '
            '25.1 °C; measurement timestamp 2026-10-07 18:23:55 UTC. Corning 3635, portrait reader profile reused from the Orange G workflow, '
            'height 10.27 mm. The profile is not independently calibrated for this plate. Original validity review remains pending.', 'small')
  paragraph('<b>Blank correction:</b> raw sample OD minus the mean of six A31 blanks at the same planned volume. '
            'Sample SD uses n-1. All original observations and negative corrected blank values are preserved. D9 is flagged without exclusion. '
            'Three wells share each sample concentration/volume; all derive from one serial dilution chain, not independent preparation replicates.', 'small')
  paragraph('<b>Uncertainty:</b> 10,000 percentile bootstrap resamples, seed 1034. Resample three sample wells within each group and six blanks '
            'within each volume; reuse the same resampled blank mean across the six concentrations sharing that volume. Fit ordinary least squares '
            'with a free intercept. Intervals are conditional on this plate and exclude between-day variation and uncertainty in concentration, dilution and actual volume.', 'small')
  paragraph('<b>Verification:</b> all 96 well mappings match the raw backend, run manifest and saved layout. The old concentration metadata is '
            'explicitly superseded by the user-reported 70 ng/µL stock. Reader blank calculations and all corrected values were reproduced. '
            'Fourteen slope/intercept/R²/RMSE results were independently checked with SciPy. Seven input hashes preserve the raw provenance. '
            'No reader or robot was connected during analysis.', 'small')
  paragraph('<b>Recommended follow-up:</b> repeat an independently prepared low-range series, for example 0.25, 0.5, 1, 2, 4 and 8 ng/µL, '
            'with matched blanks at 50, 75, 100 and 200 µL. Three samples per condition plus six blanks per volume fill 96 wells. '
            'Repeat on separate days, independently check stock concentration and delivered volume, and document bubbles or suspect wells before reading. '
            'This would test the suggested ~1 ng/µL noise limit and whether 75 µL is consistently adequate.', 'small')
  paragraph('Sources and reproducibility', 'h2')
  for i, source in enumerate(meta['sources'], 1):
    paragraph(f'[{i}] <a href="{escape(source["url"], quote=True)}" color="#007E75">{escape(source["title"])}</a>', 'small')
  paragraph('The accompanying folder contains an executed notebook, analysis and plotting scripts, original input copies with hashes, '
            'all CSV tables, bootstrap arrays and 15 PNG/SVG figure pairs. Open DNA_A260_analysis.ipynb to reproduce the numerical analysis and plots.', 'small')

  def footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor('#CDD8DC'))
    canvas.line(40, 30, width-40, 30)
    canvas.setFont('DejaVu', 8)
    canvas.setFillColor(INK)
    canvas.drawString(40, 18, 'Exp1034 DNA A260 | A30 = 70 ng/µL | exploratory, single-plate analysis')
    canvas.drawRightString(width-40, 18, str(doc.page))
    canvas.restoreState()

  out=root / 'output/pdf/DNA_A260_lab_notebook.pdf'
  out.parent.mkdir(parents=True, exist_ok=True)
  doc=SimpleDocTemplate(str(out), pagesize=landscape(A4), rightMargin=40, leftMargin=40,
                        topMargin=25, bottomMargin=40, title='Exp1034 DNA A260: volume and sensitivity', author='Exp1034 analysis')
  doc.build(story, onFirstPage=footer, onLaterPages=footer)
  css='body{max-width:1160px;margin:40px auto;padding:0 24px;font:17px/1.6 system-ui;color:#183341}h1,h2{color:#007e75}table{border-collapse:collapse;width:100%;margin:20px 0}td,th{padding:10px;text-align:left;border-bottom:1px solid #cdd8dc}th{background:#dceeeb}img{width:100%;height:auto}section{margin:45px 0}hr{border:0;border-top:1px solid #cdd8dc;margin:40px 0}.download{font-size:14px}a{color:#007e75}'
  (root / 'DNA_A260_lab_notebook.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Exp1034 DNA A260</title><style>'+css+'</style><body>'+''.join(html)+'</body></html>')
  print(out)


if __name__ == '__main__':
  build()
