# -----------------------------------------------------------------------------
# File:    rcpp/scripts/table1_simulation.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Table 1.
# Purpose: Table 1 with the Rcpp version: R, FDR and TPR of the six methods in
#          simulated data, for C = 10 and C = 5.
# Inputs:  As scripts/table1_simulation.R; rcpp/src/ (the C++ code).
# Outputs: As scripts/table1_simulation.R, same file names.
# Usage:   Rscript rcpp/scripts/table1_simulation.R [--rerun] [--runs=1000]
#          Run from the repository root. Options as scripts/table1_simulation.R,
#          which this file runs with the C++ functions in place of the R ones.
# Source:  scripts/table1_simulation.R.
# -----------------------------------------------------------------------------

source("R/setup.R")
source("rcpp/R/load.R")

run_r_script("table1_simulation")
