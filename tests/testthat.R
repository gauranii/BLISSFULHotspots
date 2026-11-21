# -----------------------------------------------------------------------------
# File:    tests/testthat.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. All figures and tables.
# Purpose: Run the unit tests of the functions in R/: every file
#          tests/testthat/test-*.R, on small synthetic inputs.
# Inputs:  tests/testthat/test-*.R; the tests read archive/ and data/domains/
#          and write only to temporary folders.
# Outputs: prints one line per test file; exit status 1 if a test fails.
# Usage:   Rscript tests/testthat.R      # from anywhere, inside `nix develop`
#          nix run .#test                # the same, with the pinned R
# -----------------------------------------------------------------------------

library(testthat)

runner <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
tests_dir <- file.path(dirname(normalizePath(runner)), "testthat")

test_dir(tests_dir, reporter = "summary", stop_on_failure = TRUE, stop_on_warning = FALSE)
