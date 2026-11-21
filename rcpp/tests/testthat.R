# -----------------------------------------------------------------------------
# File:    rcpp/tests/testthat.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. All figures and tables.
# Purpose: Run the unit tests of the Rcpp version: every file
#          rcpp/tests/testthat/test-*.R, which compare each C++ function with
#          the R function of the same name.
# Inputs:  rcpp/tests/testthat/test-*.R; the tests compile rcpp/src/, read
#          data/domains/ and write only to temporary folders.
# Outputs: prints one line per test file; exit status 1 if a test fails.
# Usage:   Rscript rcpp/tests/testthat.R   # from anywhere, inside `nix develop`
#          nix run .#test                  # with the R and the Python tests
# -----------------------------------------------------------------------------

library(testthat)

runner <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
tests_dir <- file.path(dirname(normalizePath(runner)), "testthat")

test_dir(tests_dir, reporter = "summary", stop_on_failure = TRUE, stop_on_warning = FALSE)
