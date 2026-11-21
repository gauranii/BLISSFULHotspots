# -----------------------------------------------------------------------------
# File:    rcpp/scripts/figS1_fdr_tpr_well_separated.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Figure S1.
# Purpose: Figure S1 with the Rcpp version: FDR and TPR of the six methods
#          against N when f0 and f1 are well separated.
# Inputs:  As scripts/figS1_fdr_tpr_well_separated.R; rcpp/src/ (the C++ code).
# Outputs: As scripts/figS1_fdr_tpr_well_separated.R, same file names.
# Usage:   Rscript rcpp/scripts/figS1_fdr_tpr_well_separated.R [--rerun] [--runs=1000]
#          Run from the repository root. Options as scripts/figS1_fdr_tpr_well_separated.R,
#          which this file runs with the C++ functions in place of the R ones.
# Source:  scripts/figS1_fdr_tpr_well_separated.R.
# -----------------------------------------------------------------------------

source("R/setup.R")
source("rcpp/R/load.R")

run_r_script("figS1_fdr_tpr_well_separated")
