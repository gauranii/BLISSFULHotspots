# -----------------------------------------------------------------------------
# File:    tests/testthat/helper-inputs.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Inputs shared by the unit tests (testthat loads this file before
#          the tests): temporary folders, small settings, small simulated
#          data sets and made-up domain files.
# Inputs:  None.
# Outputs: None.
# Usage:   loaded by testthat; see tests/testthat.R
# -----------------------------------------------------------------------------

#' A new empty temporary folder.
#'
#' @return Its path.
new_folder <- function() {
  folder <- tempfile("localfdr-test-")
  dir.create(folder)
  folder
}

#' Run code without what it prints.
#'
#' @param code Expression that prints time-stamped messages.
#' @return The value of `code`.
quietly <- function(code) {
  capture.output(value <- code)
  value
}

#' TRUE if a file starts as a PDF file does and is not empty of marks.
#'
#' @param file Path of the file.
#' @return TRUE or FALSE.
is_pdf <- function(file) {
  file.exists(file) && file.size(file) > 1000 &&
    identical(readBin(file, "raw", 4), charToRaw("%PDF"))
}

# Settings as PAPER of R/setup.R with short chains.
SMALL_SETTINGS <- list(alpha = 0.05, replications = 1000, sim_draws = 30, sim_burn = 10,
                       real_draws = 40, real_burn = 10, real_seed = 1234,
                       real_gamma = c(NP = 1.5, SP = 0.02, P = 0.5, GP = 0.5),
                       figS2_seed = 1439, figure_size = 3)

#' One setting in the layout of `simulation_settings()`.
#'
#' The default has the null and the non-null counts far apart: a ZIGP null
#' with mean 1.5 and non-null counts above C = 10.
#'
#' @param N,lambda,theta,pi0,C,run As the columns of `simulation_settings()`.
#' @return One-row data frame.
small_setting <- function(N = 60, lambda = 1.5, theta = 0, pi0 = 0.8, C = 10, run = "run0224") {
  data.frame(run = run, setting = 1, eta = 0.4, lambda = lambda, theta = theta, pi0 = pi0,
             f1 = "binomial", C = C, N = N)
}

#' Groups as `draw_groups()` returns them, from the null and non-null counts
#' at 0, 1, ..., K.
#'
#' @param n0,n1 Numbers of null and non-null observations at each count.
#' @param C Cutoff.
#' @return List with K, N0, N1, N2, N3, n00, n0, n1.
made_up_groups <- function(n0, n1, C = -1) {
  values <- seq_along(n0) - 1
  list(K = max(values), N0 = sum(n0), N1 = sum(n1 * (values > C)), N2 = sum(values * n0),
       N3 = sum(values * n1 * (values > C)), n00 = n0[1], n0 = n0, n1 = n1)
}

#' Write a made-up domain file in the layout of data/domains/: the number of
#' positions, then `position:count:mutations` from position 0 on.
#'
#' @param folder Folder of the file.
#' @param domain Name of the domain.
#' @param counts Mutation counts of positions 1 to N.
#' @return The path of the file.
write_domain_file <- function(folder, domain, counts) {
  file <- file.path(folder, paste0(domain, ".txt"))
  writeLines(c(length(counts), "0:0:",
               sprintf("%d:%d:m%d", seq_along(counts), counts, seq_along(counts))), file)
  file
}

# Mutation counts of two made-up domains: mostly small counts and a few
# large ones.
DOMAIN_COUNTS <- list(made1 = c(rep(0, 30), rep(1, 12), rep(2, 6), 3, 3, 14, 15, 15, 17),
                      made2 = c(rep(0, 25), rep(1, 15), rep(2, 5), 3, 12, 12, 13, 16))
