# -----------------------------------------------------------------------------
# File:    rcpp/scripts/table3_protein_domains.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Table 3.
# Purpose: Table 3 with the Rcpp version: number of positions declared hotspots
#          in the protein domains, by method.
# Inputs:  As scripts/table3_protein_domains.R; rcpp/src/ (the C++ code).
# Outputs: As scripts/table3_protein_domains.R, same file names.
# Usage:   Rscript rcpp/scripts/table3_protein_domains.R [--cores=4] [--domains=all]
#          Run from the repository root. Options as scripts/table3_protein_domains.R,
#          which this file runs with the C++ functions in place of the R ones.
# Source:  scripts/table3_protein_domains.R.
# -----------------------------------------------------------------------------

source("R/setup.R")
source("rcpp/R/load.R")

run_r_script("table3_protein_domains")
