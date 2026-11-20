# -----------------------------------------------------------------------------
# File:    checks/extract_paper_values.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Pull every number of Tables 1 to 3 out of the final manuscript
#          (main.tex) into checks/paper_values.csv.
# Usage:   python3 checks/extract_paper_values.py path/to/main.tex > checks/paper_values.csv
# Outputs: one line per cell: table, n (N, or the domain in Table 3), row
#          (method, prior, or empty), block (C10, C5, NP, SP, or empty),
#          measure (R, FDR, TPR; N, EB, NP, SP, GP in Table 3), stat (mean, sd
#          or count) and value, as printed.
# -----------------------------------------------------------------------------
import re
import sys

tex = open(sys.argv[1]).read()
MEASURES = ['R', 'FDR', 'TPR']


def tabular(label):
    """Rows of the tabular of the table with this label, as lists of cells."""
    start = tex.index('\\label{' + label + '}')
    body = tex[tex.index('\\begin{tabular}', start):tex.index('\\end{tabular}', start)]
    body = re.sub(r'\\textbf\{([^}]*)\}', r'\1', body).replace('\\hline', '')
    return [[c.strip() for c in row.split('&')] for row in body.split('\\\\')]


def number(cell):
    return re.fullmatch(r'\(?(\d+(?:\.\d+)?)\)?', cell).group(1)


print('table,n,row,block,measure,stat,value')
# Tables 1 and 2: a row of means (N, label, six numbers), then a row of
# standard deviations in parentheses.
for table, label, blocks in [(1, 'tab:T1N50', ['C10', 'C5']), (2, 'tab:P1N50', ['NP', 'SP'])]:
    n = row = ''
    for cells in tabular(label):
        if len(cells) != 8 or not re.fullmatch(r'\(?[\d.]+\)?', cells[2]):
            continue
        stat = 'sd' if cells[2].startswith('(') else 'mean'
        if stat == 'mean':
            n = cells[0] or n
            row = cells[1]
        for i, cell in enumerate(cells[2:]):
            print(','.join([str(table), n, row, blocks[i // 3], MEASURES[i % 3], stat, number(cell)]))
# Table 3: two domains per row, each with N, EB, NP, SP, GP.
for cells in tabular('tab:RejectionsB'):
    if len(cells) != 12 or not cells[1].isdigit():
        continue
    for half in (cells[:6], cells[6:]):
        for measure, cell in zip(['N', 'EB', 'NP', 'SP', 'GP'], half[1:]):
            print(','.join(['3', half[0], '', '', measure, 'count', number(cell)]))
