# -----------------------------------------------------------------------------
# File:    rcpp/scripts/table2_prior_sensitivity.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Table 2.
# Purpose: Table 2 with the Rcpp version: prior sensitivity of the nonparametric
#          and semiparametric methods.
# Inputs:  As scripts/table2_prior_sensitivity.R; rcpp/src/ (the C++ code).
# Outputs: As scripts/table2_prior_sensitivity.R, same file names.
# Usage:   Rscript rcpp/scripts/table2_prior_sensitivity.R [--rerun] [--runs=1000]
#          Run from the repository root. Options as scripts/table2_prior_sensitivity.R,
#          which this file runs with the C++ functions in place of the R ones.
# Source:  scripts/table2_prior_sensitivity.R.
# -----------------------------------------------------------------------------

source("R/setup.R")
source("rcpp/R/load.R")

run_r_script("table2_prior_sensitivity")
