# R and Python

Written by `checks/check_python.sh`: the outputs of the R scripts (`scripts/`)
and of the Python scripts (`python/scripts/`) compared by code
(`checks/check_python.py`). Every csv file under `results/` is compared cell
by cell, without tolerance: the same header, the same rows, every cell the
same text or the same number. Every figure is compared by the marks read from
the two PDF files (`checks/figure_compare.py`): the same number of bars, lines
or tiles in the same order and with the same fills; bar heights, bar positions
and line vertices agree to 0.02 pt of the R page after one straight-line map
per axis, since R and matplotlib lay their pages out differently and R prints
its coordinates with two decimals. "Values" counts the compared cells and
mark coordinates. Times are wall-clock seconds on an 8 core Apple M2 laptop
shared with other jobs; the Python times include loading its compiled (numba)
functions.

The tests in `python/tests/` (also run by `nix flake check`) compare with R,
bit for bit: the random numbers of `python/localfdr/rng.py` (uniform, normal,
exponential, gamma, beta, binomial, Poisson, geometric, `sample.int`) for
several seeds and over every branch of the algorithms; `lgamma`, the logistic
functions and `round`; the Cholesky factor, the eigenvalues and the
covariances behind the adapted proposal; the ZIGP probabilities and draws;
the empirical Bayes fit of five simulated data sets and of the six domains;
and every stored draw of seeded Gibbs chains for the three versions of the
sampler, the four alternatives and three data sets.

## With the defaults of the scripts

Both implementations run with the defaults of their scripts, as
`nix flake check` runs them: Tables 1 and 2 and Figures S1 and S3 from the
archived runs, Table 3 from the data (18 seeded chains of 11000 iterations),
Figure 1 from the data and Figure S2 from seeded simulated data.

7 of 7 outputs identical in every value (2 cores).

| Output | From | Files | Values | Not identical | Largest difference | Verdict | R (s) | Python (s) |
|---|---|---|---|---|---|---|---|---|
| Table 1 | archived runs | 1 | 420 | 0 | 0 | identical | 0.2 | 0.4 |
| Table 2 | archived runs | 1 | 375 | 0 | 0 | identical | 0.2 | 0.5 |
| Table 3 | computed (seeded chains) | 1 | 36 | 0 | 0 | identical | 55.3 | 10.9 |
| Figure 1 | computed | 7 | 1524 | 0 | 0 | identical | 0.2 | 1.2 |
| Figure S1 | archived runs | 5 | 1848 | 0 | 0 | identical | 1.0 | 1.0 |
| Figure S2 | computed (seeded data) | 10 | 3245 | 0 | 0 | identical | 0.2 | 1.6 |
| Figure S3 | archived runs | 2 | 1480 | 0 | 0 | identical | 0.6 | 1.3 |

## Recomputed by simulation

Both implementations run Tables 1 and 2 with `--rerun` for the first
replications of every setting and method, with the seeds of the full runs
(`bash checks/check_python.sh rerun <work dir>`). Compared are the tables and
the stored file of every replication (rejections, false and true rejections
of each data set and method).

2 of 2 outputs identical in every value (2 cores).

| Output | From | Files | Values | Not identical | Largest difference | Verdict | R (s) | Python (s) |
|---|---|---|---|---|---|---|---|---|
| Table 1 | 2 replications | 81 | 940 | 0 | 0 | identical | 100.1 | 19.2 |
| Table 2 | 2 replications | 101 | 875 | 0 | 0 | identical | 126.5 | 27.4 |

