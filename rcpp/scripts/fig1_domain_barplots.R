# -----------------------------------------------------------------------------
# File:    rcpp/scripts/fig1_domain_barplots.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Figure 1.
# Purpose: Figure 1 with the Rcpp version loaded: bar plots of six protein
#          domains. The figure runs no sampler, so it is the R figure.
# Inputs:  As scripts/fig1_domain_barplots.R; rcpp/src/ (the C++ code).
# Outputs: As scripts/fig1_domain_barplots.R, same file names.
# Usage:   Rscript rcpp/scripts/fig1_domain_barplots.R
#          Run from the repository root. Options as scripts/fig1_domain_barplots.R,
#          which this file runs with the C++ functions in place of the R ones.
# Source:  scripts/fig1_domain_barplots.R.
# -----------------------------------------------------------------------------

source("R/setup.R")
source("rcpp/R/load.R")

run_r_script("fig1_domain_barplots")
