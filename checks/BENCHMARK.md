# R and Rcpp run times

Written by `checks/write_benchmark.R` from `bash checks/benchmark.sh <work dir>`
(see Rerun). Times are wall-clock seconds on an 8 core Apple M2 laptop (16 GB), with no other heavy jobs,
OpenBLAS on one thread, the C++ files compiled beforehand. The chain is timed
inside R; the scripts are timed whole, start of R included. Every run of R and
of Rcpp gives the same output files, byte for byte (last column).

## Same seeded work, 5 runs each

Median (range) over the runs; speedup = median R / median Rcpp.

| Work | Cores | Runs | R (s) | Rcpp (s) | Speedup | Outputs identical |
|---|---|---|---|---|---|---|
| One Gibbs chain (cd00063, NP, 11000 iterations) | 1 | 5 | 5.81 (5.45 to 6.01) | 0.35 (0.34 to 0.36) | 16.6x | 5 of 5 |
| Table 3 (18 chains of 11000 iterations) | 1 | 5 | 80.4 (78.8 to 85.1) | 6.70 (6.63 to 6.83) | 12.0x | 5 of 5 |
| Table 3 (18 chains of 11000 iterations) | 2 | 5 | 56.9 (54.9 to 69.5) | 4.50 (4.46 to 5.24) | 12.7x | 5 of 5 |
| Table 1, one replication of every setting and method | 1 | 5 | 75.2 (73.7 to 78.7) | 7.97 (7.90 to 8.05) | 9.4x | 5 of 5 |
| Table 1, one replication of every setting and method | 2 | 5 | 52.9 (49.9 to 66.6) | 5.27 (5.10 to 6.46) | 10.0x | 5 of 5 |

## The full simulation runs, estimated

One replication of every setting and method of each simulation study, run
once on 1 core, in the order below (replications shared by two studies are
run once, as in a full run: Table 2 reuses the uniform-prior replications of
Table 1). The full runs have 1000 replications; the estimate is the 1-core
time of one replication times 1000, divided by 8 cores, as if the 8 cores
shared the work evenly.

| Study | R, 1 replication (s) | Rcpp (s) | Speedup | R, full run, 8 cores | Rcpp, full run, 8 cores |
|---|---|---|---|---|---|
| Table 1 | 74.0 | 8.20 | 9.0x | 2.6 h | 17 min |
| Table 2 | 87.2 | 5.89 | 14.8x | 3.0 h | 12 min |
| Figure S1 | 186.9 | 19.4 | 9.6x | 6.5 h | 40 min |
| Figure S3 | 489.3 | 44.9 | 10.9x | 17.0 h | 1.6 h |
| All four | 837.4 | 78.4 | 10.7x | 29.1 h | 2.7 h |

Outputs of the two versions identical: yes (392 files).

## Rerun

```sh
nix develop
bash checks/benchmark.sh /tmp/benchmark          # RUNS=5 by default; about 42 min
Rscript checks/write_benchmark.R /tmp/benchmark "<machine>" > checks/BENCHMARK.md
```
