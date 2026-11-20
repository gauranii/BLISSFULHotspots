# Check results

Written by `checks/run_checks.sh`: the outputs of the scripts in `scripts/`
(run with their defaults) compared with the paper by code.

## Tables 1 to 3

Every cell is compared with the number printed in the paper
(`checks/paper_values.csv`).

| Table | Cells equal | Cells | Check |
|---|---|---|---|
| Table 1 (simulation, six methods) | 360 | 360 | OK |
| Table 2 (prior sensitivity) | 300 | 300 | OK |
| Table 3 (protein domains) | 30 | 30 | OK |

### By row

Cells equal / cells of each row group: Tables 1 and 2 by N (means and
standard deviations of R, FDR and TPR), Table 3 by domain (N, EB, NP, SP, GP).

| Table | Row group | Cells equal | Cells | Check |
|---|---|---|---|---|
| 1 | N = 50 | 72 | 72 | OK |
| 1 | N = 100 | 72 | 72 | OK |
| 1 | N = 200 | 72 | 72 | OK |
| 1 | N = 500 | 72 | 72 | OK |
| 1 | N = 1000 | 72 | 72 | OK |
| 2 | N = 50 | 60 | 60 | OK |
| 2 | N = 100 | 60 | 60 | OK |
| 2 | N = 200 | 60 | 60 | OK |
| 2 | N = 500 | 60 | 60 | OK |
| 2 | N = 1000 | 60 | 60 | OK |
| 3 | cd00054 | 5 | 5 | OK |
| 3 | cd00063 | 5 | 5 | OK |
| 3 | cd00083 | 5 | 5 | OK |
| 3 | cd00590 | 5 | 5 | OK |
| 3 | smart00179 | 5 | 5 | OK |
| 3 | smart00209 | 5 | 5 | OK |

### Table 3, cell by cell

Each cell reads paper / computed.

| Domain | N | EB | NP | SP | GP |
|---|---|---|---|---|---|
| cd00054 | 90 / 90 | 45 / 45 | 43 / 43 | 42 / 42 | 40 / 40 |
| cd00063 | 193 / 193 | 100 / 100 | 95 / 95 | 89 / 89 | 89 / 89 |
| cd00083 | 93 / 93 | 7 / 7 | 1 / 1 | 0 / 0 | 1 / 1 |
| cd00590 | 175 / 175 | 74 / 74 | 69 / 69 | 19 / 19 | 66 / 66 |
| smart00179 | 182 / 182 | 44 / 44 | 43 / 43 | 38 / 38 | 38 / 38 |
| smart00209 | 188 / 188 | 54 / 54 | 48 / 48 | 14 / 14 | 32 / 32 |

## Figures 1, S1, S2 and S3

Figures 1, S1 and S2 are compared with the published files (`archive/figures/`) by
the marks in the two PDF files, Figure S1 also by its plotted values, and Figure S3
by the cells of its table. 24 of 24 checks OK.

| Figure | Check | Equal | Of | Detail | Result |
|---|---|---|---|---|---|
| Fig 1A | rectangles equal to the published file | 53 | 53 | largest difference 0.00 pt | OK |
| Fig 1B | rectangles equal to the published file | 22 | 22 | largest difference 0.00 pt | OK |
| Fig 1C | rectangles equal to the published file | 29 | 29 | largest difference 0.00 pt | OK |
| Fig 1D | rectangles equal to the published file | 48 | 48 | largest difference 0.00 pt | OK |
| Fig 1E | rectangles equal to the published file | 46 | 46 | largest difference 0.00 pt | OK |
| Fig 1F | rectangles equal to the published file | 56 | 56 | largest difference 0.00 pt | OK |
| Fig S1A | lines equal to the published file | 6 | 6 | largest difference 0.00 pt | OK |
| Fig S1B | lines equal to the published file | 6 | 6 | largest difference 0.00 pt | OK |
| Fig S1C | lines equal to the published file | 6 | 6 | largest difference 0.00 pt | OK |
| Fig S1D | lines equal to the published file | 6 | 6 | largest difference 0.00 pt | OK |
| Fig S2A | rectangles equal to the published file | 53 | 53 | largest difference 0.00 pt | OK |
| Fig S2B | rectangles equal to the published file | 69 | 69 | largest difference 0.00 pt | OK |
| Fig S2C | rectangles equal to the published file | 69 | 69 | largest difference 0.00 pt | OK |
| Fig S2D | rectangles equal to the published file | 55 | 55 | largest difference 0.00 pt | OK |
| Fig S2E | rectangles equal to the published file | 73 | 73 | largest difference 0.00 pt | OK |
| Fig S2F | rectangles equal to the published file | 73 | 73 | largest difference 0.00 pt | OK |
| Fig S2G | rectangles equal to the published file | 57 | 57 | largest difference 0.00 pt | OK |
| Fig S2H | rectangles equal to the published file | 75 | 75 | largest difference 0.00 pt | OK |
| Fig S2I | rectangles equal to the published file | 75 | 75 | largest difference 0.00 pt | OK |
| Fig S1 | FDR, setting 16: plotted values equal to the plot input | 66 | 66 |  | OK |
| Fig S1 | TPR, setting 16: plotted values equal to the plot input | 66 | 66 |  | OK |
| Fig S1 | FDR, setting 18: plotted values equal to the plot input | 66 | 66 |  | OK |
| Fig S1 | TPR, setting 18: plotted values equal to the plot input | 66 | 66 |  | OK |
| Fig S3 | cells equal to the sheet of the published table | 360 | 360 |  | OK |
