# R and Rcpp

Written by `checks/check_rcpp.sh`: the outputs of the R scripts (`scripts/`)
and of the Rcpp scripts (`rcpp/scripts/`) compared by code
(`checks/check_python.py`, the comparison also used for the Python version).
Every csv file under `results/` is compared cell by cell, without tolerance:
the same header, the same rows, every cell the same text or the same number.
Every figure is compared by the marks read from the two PDF files
(`checks/figure_compare.py`): the same number of bars, lines or tiles in the
same order, with the same fills and positions. "Values" counts the compared
cells and mark coordinates. Times are wall-clock seconds on an 8 core Apple
M2 laptop shared with other jobs.
Compiling the C++ files took 8 s, once; the Rcpp times are those of
the scripts with the compiled files in place.

The Rcpp scripts run the R scripts with the C++ functions of `rcpp/src/` in
place of the R functions of the same names: the distributions, the whole
Gibbs sampler and the simulated data. The figures themselves, the empirical
Bayes methods and the tables are made by the same R code in both versions.

The tests in `rcpp/tests/` (also run by `nix flake check`) call every C++
function and its R function with the same arguments under the same seed and
expect identical values and the same state of the random number generator
afterwards: the ZIGP probabilities and draws, the Dirichlet and multinomial
draws, every step of the sampler, the adapted proposal, the simulated data,
and seeded Gibbs chains for the three versions of the sampler, the four
alternatives, both switches and three data sets.

## With the defaults of the scripts

Both versions run with the defaults of their scripts, as `nix flake check`
runs them: Tables 1 and 2 and Figures S1 and S3 from the archived runs,
Table 3 from the data (18 seeded chains of 11000 iterations), Figure 1 from
the data and Figure S2 from seeded simulated data.

7 of 7 outputs identical in every value (2 cores).

| Output | From | Files | Values | Not identical | Largest difference | Verdict | R (s) | Rcpp (s) |
|---|---|---|---|---|---|---|---|---|
| Table 1 | archived runs | 1 | 420 | 0 | 0 | identical | 0.2 | 0.3 |
| Table 2 | archived runs | 1 | 375 | 0 | 0 | identical | 0.2 | 0.3 |
| Table 3 | computed (seeded chains) | 1 | 36 | 0 | 0 | identical | 54.4 | 4.5 |
| Figure 1 | computed | 7 | 1524 | 0 | 0 | identical | 0.2 | 0.3 |
| Figure S1 | archived runs | 5 | 1848 | 0 | 0 | identical | 0.9 | 1.0 |
| Figure S2 | computed (seeded data) | 10 | 3245 | 0 | 0 | identical | 0.2 | 0.3 |
| Figure S3 | archived runs | 2 | 1480 | 0 | 0 | identical | 0.6 | 0.7 |

## Recomputed by simulation

Both versions run Tables 1 and 2 with `--rerun` for the first replications of
every setting and method, with the seeds of the full runs
(`bash checks/check_rcpp.sh rerun <work dir>`). Compared are the tables and
the stored file of every replication (rejections, false and true rejections
of each data set and method).

2 of 2 outputs identical in every value (2 cores).

| Output | From | Files | Values | Not identical | Largest difference | Verdict | R (s) | Rcpp (s) |
|---|---|---|---|---|---|---|---|---|
| Table 1 | 2 replications | 81 | 940 | 0 | 0 | identical | 100.2 | 9.9 |
| Table 2 | 2 replications | 101 | 875 | 0 | 0 | identical | 119.2 | 6.6 |

