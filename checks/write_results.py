# -----------------------------------------------------------------------------
# File:    checks/write_results.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Turn the outputs of check_tables.py and check_figures.py into
#          checks/RESULTS.md.
# Usage:   python3 checks/write_results.py <checks output dir> > checks/RESULTS.md
# -----------------------------------------------------------------------------
import csv
import os
import sys

out = sys.argv[1]
tables = list(csv.DictReader(open(os.path.join(out, 'tables.csv'))))
figures = list(csv.DictReader(open(os.path.join(out, 'figures.csv'))))
NAMES = {'1': 'Table 1 (simulation, six methods)', '2': 'Table 2 (prior sensitivity)',
         '3': 'Table 3 (protein domains)'}


def mark(ok):
    return 'OK' if ok else 'FAIL'


def count(cells):
    return sum(c['ok'] == 'True' for c in cells), len(cells)


md = ['# Check results', '',
      'Written by `checks/run_checks.sh`: the outputs of the scripts in `scripts/`',
      '(run with their defaults) compared with the paper by code.', '',
      '## Tables 1 to 3', '',
      'Every cell is compared with the number printed in the paper',
      '(`checks/paper_values.csv`).', '',
      '| Table | Cells equal | Cells | Check |', '|---|---|---|---|']
for n in '123':
    equal, total = count([c for c in tables if c['table'] == n])
    md.append(f'| {NAMES[n]} | {equal} | {total} | {mark(equal == total)} |')
md += ['', '### By row', '',
       'Cells equal / cells of each row group: Tables 1 and 2 by N (means and',
       'standard deviations of R, FDR and TPR), Table 3 by domain (N, EB, NP, SP, GP).', '',
       '| Table | Row group | Cells equal | Cells | Check |', '|---|---|---|---|---|']
for n in '123':
    groups = []
    for c in tables:
        if c['table'] == n and c['n'] not in groups:
            groups.append(c['n'])
    for g in groups:
        equal, total = count([c for c in tables if c['table'] == n and c['n'] == g])
        label = g if n == '3' else f'N = {g}'
        md.append(f'| {n} | {label} | {equal} | {total} | {mark(equal == total)} |')
md += ['', '### Table 3, cell by cell', '', 'Each cell reads paper / computed.', '',
       '| Domain | N | EB | NP | SP | GP |', '|---|---|---|---|---|---|']
domains = []
for c in tables:
    if c['table'] == '3' and c['n'] not in domains:
        domains.append(c['n'])
for d in domains:
    cells = {c['measure']: f"{c['value']} / {c['computed']}" for c in tables
             if c['table'] == '3' and c['n'] == d}
    md.append(f"| {d} | " + ' | '.join(cells[m] for m in ['N', 'EB', 'NP', 'SP', 'GP']) + ' |')

equal, total = count(figures)
md += ['', '## Figures 1, S1, S2 and S3', '',
       'Figures 1, S1 and S2 are compared with the published files (`archive/figures/`) by',
       'the marks in the two PDF files, Figure S1 also by its plotted values, and Figure S3',
       f'by the cells of its table. {equal} of {total} checks OK.', '',
       '| Figure | Check | Equal | Of | Detail | Result |', '|---|---|---|---|---|---|']
md += [f"| {r['figure']} | {r['check']} | {r['equal']} | {r['total']} | {r['detail']} | "
       f"{mark(r['ok'] == 'True')} |" for r in figures]
print('\n'.join(md))
