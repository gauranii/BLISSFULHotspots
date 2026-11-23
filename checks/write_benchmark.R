# -----------------------------------------------------------------------------
# File:    checks/write_benchmark.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Write checks/BENCHMARK.md from the times and comparisons of
#          checks/benchmark.sh: median and range of each cell in R and Rcpp,
#          the speedup, and the estimated time of the full simulation runs.
# Inputs:  <work dir>/times.csv, <work dir>/identical.csv.
# Outputs: Markdown on standard output.
# Usage:   Rscript checks/write_benchmark.R <work dir> "<machine>" > checks/BENCHMARK.md
# Source:  checks/benchmark.sh.
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
work <- args[1]
machine <- args[2]
times <- read.csv(file.path(work, "times.csv"))
same <- read.csv(file.path(work, "identical.csv"))

fmt <- function(s) ifelse(s < 10, sprintf("%.2f", s), sprintf("%.1f", s))
hours <- function(s) if (s >= 3600) sprintf("%.1f h", s / 3600) else sprintf("%.0f min", s / 60)

cells <- data.frame(
  cell = c("chain", "table3", "table3", "table1", "table1"),
  cores = c(1, 1, 2, 1, 2),
  label = c("One Gibbs chain (cd00063, NP, 11000 iterations)",
            "Table 3 (18 chains of 11000 iterations)", "Table 3 (18 chains of 11000 iterations)",
            "Table 1, one replication of every setting and method",
            "Table 1, one replication of every setting and method"))
rows <- lapply(seq_len(nrow(cells)), function(i) {
  pick <- function(impl) times$seconds[times$cell == cells$cell[i] & times$cores == cells$cores[i] &
                                         times$impl == impl]
  r <- pick("r"); cpp <- pick("rcpp")
  ok <- same[same$cell == cells$cell[i] & same$cores == cells$cores[i], ]
  sprintf("| %s | %d | %d | %s (%s to %s) | %s (%s to %s) | %.1fx | %d of %d |",
          cells$label[i], cells$cores[i], length(r),
          fmt(median(r)), fmt(min(r)), fmt(max(r)),
          fmt(median(cpp)), fmt(min(cpp)), fmt(max(cpp)),
          median(r) / median(cpp), sum(ok$identical == "yes"), nrow(ok))
})

studies <- c(table1_simulation = "Table 1", table2_prior_sensitivity = "Table 2",
             figS1_fdr_tpr_well_separated = "Figure S1", figS3_fdr_table = "Figure S3")
study_rows <- lapply(names(studies), function(s) {
  pick <- function(impl) times$seconds[times$cell == paste0("study-", s) & times$impl == impl]
  r <- pick("r"); cpp <- pick("rcpp")
  sprintf("| %s | %s | %s | %.1fx | %s | %s |", studies[[s]], fmt(r), fmt(cpp), r / cpp,
          hours(r * 1000 / 8), hours(cpp * 1000 / 8))
})
total <- function(impl) sum(times$seconds[startsWith(times$cell, "study-") & times$impl == impl])
tr <- total("r"); tc <- total("rcpp")
study_ok <- same[same$cell == "study", ]

cat(sprintf("# R and Rcpp run times

Written by `checks/write_benchmark.R` from `bash checks/benchmark.sh <work dir>`
(see Rerun). Times are wall-clock seconds on %s, with no other heavy jobs,
OpenBLAS on one thread, the C++ files compiled beforehand. The chain is timed
inside R; the scripts are timed whole, start of R included. Every run of R and
of Rcpp gives the same output files, byte for byte (last column).

## Same seeded work, %d runs each

Median (range) over the runs; speedup = median R / median Rcpp.

| Work | Cores | Runs | R (s) | Rcpp (s) | Speedup | Outputs identical |
|---|---|---|---|---|---|---|
%s

## The full simulation runs, estimated

One replication of every setting and method of each simulation study, run
once on 1 core, in the order below (replications shared by two studies are
run once, as in a full run: Table 2 reuses the uniform-prior replications of
Table 1). The full runs have 1000 replications; the estimate is the 1-core
time of one replication times 1000, divided by 8 cores, as if the 8 cores
shared the work evenly.

| Study | R, 1 replication (s) | Rcpp (s) | Speedup | R, full run, 8 cores | Rcpp, full run, 8 cores |
|---|---|---|---|---|---|
%s
| All four | %s | %s | %.1fx | %s | %s |

Outputs of the two versions identical: %s (%d files).

## Rerun

```sh
nix develop
bash checks/benchmark.sh /tmp/benchmark          # RUNS=5 by default; about %s
Rscript checks/write_benchmark.R /tmp/benchmark \"<machine>\" > checks/BENCHMARK.md
```
",
  machine, max(times$run), paste(unlist(rows), collapse = "\n"),
  paste(unlist(study_rows), collapse = "\n"),
  fmt(tr), fmt(tc), tr / tc, hours(tr * 1000 / 8), hours(tc * 1000 / 8),
  ifelse(all(study_ok$identical == "yes"), "yes", "NO"), sum(study_ok$files),
  hours(sum(times$seconds) + 60)))
