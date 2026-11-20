# -----------------------------------------------------------------------------
# File:    R/distributions.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Tables 1 to 3, Figures S1 to S3.
# Purpose: The zero-inflated generalized Poisson (ZIGP) null distribution, the
#          shifted Poisson and generalized Poisson alternatives, and the
#          Dirichlet and multinomial samplers of the Gibbs sampler.
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R")
# Source:  simFDR_fcns.R (2021 simulation runs), fcns1_dist.R (2018 analysis of
#          the protein domains): dzigp, rzigp, rdirichlet, shiftf, rmultinomial.
# -----------------------------------------------------------------------------

#' Probability mass function of the zero-inflated generalized Poisson.
#'
#' P(X = 0) = eta + (1 - eta) exp(-lambda) and, for x > 0,
#' P(X = x) = (1 - eta) lambda (lambda + theta x)^(x - 1) exp(-lambda - theta x) / x!.
#' theta = 0 gives the zero-inflated Poisson.
#'
#' @param x Counts (non-negative).
#' @param eta Zero-inflation probability, in [0, 1].
#' @param lambda Rate, positive.
#' @param theta Dispersion, in [0, 1].
#' @param check If TRUE stop when a parameter is outside its range.
#' @return Vector of probabilities, as long as `x`.
dzigp <- function(x, eta, lambda, theta, check = TRUE) {
  if (check) {
    stopifnot(eta >= 0 && eta <= 1)
    stopifnot(lambda > 0)
    stopifnot(theta >= 0 && theta <= 1)
  }
  gp <- exp(log(lambda) + (x - 1) * log(lambda + theta * x) - lgamma(x + 1) - (lambda + theta * x))
  density <- numeric(length(x))
  zero <- which(x == 0)
  positive <- which(x != 0)
  density[zero] <- eta + (1 - eta) * gp[zero]
  density[positive] <- (1 - eta) * gp[positive]
  density
}

#' One draw from the zero-inflated generalized Poisson, by inversion.
#'
#' Uses one uniform number. The probabilities are summed from zero upwards
#' until they pass the uniform number; if a probability overflows it is
#' continued on the log scale from the previous one.
#'
#' @param eta,lambda,theta Parameters, as in `dzigp()`.
#' @param check If TRUE stop when a parameter is outside its range.
#' @return One count.
rzigp <- function(eta, lambda, theta, check = TRUE) {
  if (check) {
    stopifnot(eta >= 0 && eta <= 1)
    stopifnot(lambda > 0)
    stopifnot(theta >= 0 && theta <= 1)
  }
  block <- 1000
  cdf <- numeric(block)
  unif <- runif(1)
  cdf[1] <- eta + (1 - eta) * exp(-lambda)

  if (unif > 0) {
    recursive <- FALSE
    j <- 1
    while (cdf[j] < unif) {
      if (!recursive) {
        p <- exp(log(lambda) + (j - 1) * log(lambda + theta * j) - lgamma(j + 1) -
                   (lambda + theta * j) + log(1 - eta))
      }
      if (p == Inf) {
        recursive <- TRUE
        log_p_previous <- log(lambda) + (j - 2) * log(lambda + theta * (j - 1)) -
          lgamma(j) - (lambda + theta * (j - 1)) + log(1 - eta)
      }
      if (recursive) {
        log_p <- log_p_previous + log(lambda + theta * (j - 1)) - log(j) - theta +
          (j - 1) * log(lambda + theta * j) - (j - 1) * log(lambda + theta * (j - 1))
        log_p_previous <- log_p
        p <- exp(log_p)
      }
      if (j %% block == 0) cdf <- c(cdf, unif)
      cdf[j + 1] <- cdf[j] + p
      j <- j + 1
    }
  }

  j <- 1
  while (unif > cdf[j]) j <- j + 1
  j - 1
}

#' One draw from a Dirichlet distribution, from independent gamma draws.
#'
#' @param beta Concentration parameters.
#' @param replace_missing If TRUE a missing concentration is set to 1e-5.
#' @return Probabilities that sum to one, as long as `beta`.
rdirichlet <- function(beta, replace_missing = FALSE) {
  if (replace_missing) beta[is.na(beta)] <- 0.00001
  draw <- rgamma(length(beta), beta, 1)
  total <- 0
  for (value in draw) total <- total + value  # summed in order, in double precision
  draw / total
}

#' One draw from a multinomial distribution, as a chain of binomial draws.
#'
#' @param size Number of trials.
#' @param prob Probabilities of the categories (need not sum to one).
#' @return Counts per category.
rmultinomial <- function(size, prob) {
  K <- length(prob)
  draw <- numeric(K)
  remaining <- sum(prob)
  for (i in seq_along(prob)) {
    if (prob[i] != 0) {
      p <- prob[i] / remaining
      draw[i] <- if (p < 1) rbinom(1, size, p) else size
      size <- size - draw[i]
    }
    if (size <= 0) return(draw)
    remaining <- remaining - prob[i]
  }
  draw[K - 1] <- size
  draw
}

#' Parametric alternative density f1: a Poisson or generalized Poisson shifted
#' to start right above the cutoff C.
#'
#' @param model "P" (Poisson with mean `xi`) or "GP" (generalized Poisson with
#'   rate `delta` and dispersion `nu`).
#' @param x Counts.
#' @param xi,delta,nu Parameters of f1.
#' @param C Cutoff: f1 is zero for counts up to C.
#' @return Vector of f1(x).
shifted_density <- function(model, x, xi, delta, nu, C) {
  f1 <- numeric(length(x))
  shift <- x - C - 1
  above <- which(x > C)
  if (model == "P") {
    f1[above] <- exp(shift[above] * log(xi) - lgamma(shift[above] + 1) - xi)
  } else if (model == "GP") {
    f1[above] <- exp(log(delta) - lgamma(shift[above] + 1) +
                       (shift[above] - 1) * log(delta + nu * shift[above]) -
                       (delta + nu * shift[above]))
  }
  f1
}
