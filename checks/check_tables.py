# -----------------------------------------------------------------------------
# File:    checks/check_tables.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Compare every cell of Tables 1 to 3 written by the table scripts
#          with the number printed in the paper.
# Usage:   python3 checks/check_tables.py <results dir> <checks output dir>
# Inputs:  checks/paper_values.csv, <results dir>/Table<N>_*.csv.
# Outputs: <checks output dir>/tables.csv; prints one line per table; exit
#          status 1 unless every cell is equal.
# -----------------------------------------------------------------------------
import csv
import os
import sys

results, out = sys.argv[1:3]
here = os.path.dirname(os.path.abspath(__file__))
FILES = {'1': 'Table1_simulation', '2': 'Table2_prior-sensitivity', '3': 'Table3_protein-domains'}


def read(path):
    return list(csv.DictReader(open(path)))


def computed_cell(cell):
    """The cell of a built table that a cell of the paper corresponds to."""
    table = tables[cell['table']]
    if cell['table'] == '3':
        row = [r for r in table if r['domain'] == cell['n']]
        column = cell['measure']
    else:
        label = 'method' if cell['table'] == '1' else 'prior'
        row = [r for r in table if r['N'] == cell['n'] and r[label] == cell['row']]
        column = f"{cell['block']}_{cell['measure']}" + ('_sd' if cell['stat'] == 'sd' else '')
    assert len(row) == 1
    return row[0][column]


tables = {n: read(os.path.join(results, name + '.csv')) for n, name in FILES.items()}
rows = []
for cell in read(os.path.join(here, 'paper_values.csv')):
    computed = computed_cell(cell)
    rows.append(dict(cell, computed=computed,
                     ok=abs(float(computed) - float(cell['value'])) < 1e-9))
del cell

os.makedirs(out, exist_ok=True)
with open(os.path.join(out, 'tables.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
for n in FILES:
    cells = [r for r in rows if r['table'] == n]
    print(f"Table {n}: {sum(r['ok'] for r in cells)}/{len(cells)} cells equal to the paper")
sys.exit(0 if all(r['ok'] for r in rows) else 1)
