# -----------------------------------------------------------------------------
# File:    checks/check_figures.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Compare every figure built by the scripts with the published one.
#          Figures 1 and S2: every rectangle of each PDF with the published
#          file (position and size, to 0.01 pt). Figure S1: the plotted values
#          with the table the published plots were drawn from, and the six
#          lines of each PDF with those of the published file (every vertex,
#          to 0.01 pt). Figure S3: every cell of the table with the sheet the
#          published table shows.
# Usage:   python3 checks/check_figures.py <published dir> <figures dir>
#          <results dir> <checks output dir>
# Inputs:  <published dir>/Fig*.pdf (archive/figures), <figures dir>/Fig*.pdf,
#          <results dir>/FigS1_values.csv, FigS3_fdr_table.csv,
#          archive/figS1_plot_input.csv, archive/figS3_full_format.csv.
# Outputs: <checks output dir>/figures.csv; prints one line per check; exit
#          status 1 unless every check is OK.
# -----------------------------------------------------------------------------
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pdf_marks import TOLERANCE, polylines, rectangles, same_marks  # noqa: E402

published, drawn, results, out = sys.argv[1:5]
archive = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'archive')
rows = []


def report(figure, check, equal, total, detail=''):
    ok = equal == total and total > 0
    rows.append({'figure': figure, 'check': check, 'equal': equal, 'total': total,
                 'detail': detail, 'ok': ok})
    print(f"{figure:10} {check:52} {equal:>4}/{total:<4} {detail:22} {'OK' if ok else 'FAIL'}")


def read(path):
    return list(csv.DictReader(open(path)))


def marks(figure, name, a, b, check):
    na, nb, diff = same_marks(a, b)
    same = diff is not None and diff <= TOLERANCE
    detail = 'counts differ' if diff is None else f'largest difference {diff:.2f} pt'
    report(figure, check, na if same else 0, nb, detail)


# A drawn figure has the name of its published file.
names = sorted(n for n in os.listdir(published) if n.endswith('.pdf'))
for name in names:
    figure = 'Fig ' + name[3:name.index('_')]
    a, b = os.path.join(drawn, name), os.path.join(published, name)
    if name.startswith(('Fig1', 'FigS2')):
        marks(figure, name, rectangles(a), rectangles(b), 'rectangles equal to the published file')
    elif name.startswith('FigS1'):
        marks(figure, name, polylines(a, 11), polylines(b, 11), 'lines equal to the published file')

# Figure S1: plotted values against the tables the published plots were drawn from.
plotted = {(r['setting'], r['measure'], r['N'], r['method']): float(r['value'])
           for r in read(os.path.join(results, 'FigS1_values.csv'))}
source = read(os.path.join(archive, 'figS1_plot_input.csv'))
methods = [c for c in source[0] if c not in ('setting', 'measure', 'C', 'N')]
for setting, measure in sorted({(r['setting'], r['measure']) for r in source}):
    cells = [(plotted.get((setting, measure, r['N'], m)), float(r[m]))
             for r in source if (r['setting'], r['measure']) == (setting, measure) for m in methods]
    report('Fig S1', f'{measure}, setting {setting}: plotted values equal to the plot input',
           sum(a is not None and abs(a - b) < 1e-9 for a, b in cells), len(cells))

# Figure S3: every cell against the sheet "Full Format".
table = {(r['N'], r['lambda']): r for r in read(os.path.join(results, 'FigS3_fdr_table.csv'))}
sheet = read(os.path.join(archive, 'figS3_full_format.csv'))
equal = sum(abs(float(table[(r['N'], r['lambda'])][f"{r['method']}_theta{r['theta']}"]) -
                float(r['FDR'])) < 1e-9 for r in sheet)
report('Fig S3', 'cells equal to the sheet of the published table', equal, len(sheet))

os.makedirs(out, exist_ok=True)
with open(os.path.join(out, 'figures.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
sys.exit(0 if all(r['ok'] for r in rows) else 1)
