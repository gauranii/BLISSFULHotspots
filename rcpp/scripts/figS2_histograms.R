# -----------------------------------------------------------------------------
# File:    rcpp/scripts/figS2_histograms.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Figure S2.
# Purpose: Figure S2 with the Rcpp version: histograms of nine simulated data
#          sets, drawn by rmixture() in C++.
# Inputs:  As scripts/figS2_histograms.R; rcpp/src/ (the C++ code).
# Outputs: As scripts/figS2_histograms.R, same file names.
# Usage:   Rscript rcpp/scripts/figS2_histograms.R
#          Run from the repository root. Options as scripts/figS2_histograms.R,
#          which this file runs with the C++ functions in place of the R ones.
# Source:  scripts/figS2_histograms.R.
# -----------------------------------------------------------------------------

source("R/setup.R")
source("rcpp/R/load.R")

run_r_script("figS2_histograms")
