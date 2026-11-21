# -----------------------------------------------------------------------------
# File:    rcpp/scripts/figS3_fdr_table.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Figure S3.
# Purpose: Figure S3 with the Rcpp version: table of the FDR of the six methods,
#          shaded by FDR.
# Inputs:  As scripts/figS3_fdr_table.R; rcpp/src/ (the C++ code).
# Outputs: As scripts/figS3_fdr_table.R, same file names.
# Usage:   Rscript rcpp/scripts/figS3_fdr_table.R [--rerun] [--runs=1000]
#          Run from the repository root. Options as scripts/figS3_fdr_table.R,
#          which this file runs with the C++ functions in place of the R ones.
# Source:  scripts/figS3_fdr_table.R.
# -----------------------------------------------------------------------------

source("R/setup.R")
source("rcpp/R/load.R")

run_r_script("figS3_fdr_table")
