# -----------------------------------------------------------------------------
# File:    R/eb_fdr.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Tables 1 and 3, Figures S1 and S3.
# Purpose: The empirical Bayes competitors: the ZIGP null fitted by an EM
#          algorithm to the counts up to a cutoff, the choice of the cutoff,
#          and the three procedures built on the fit: the two-stage and the
#          one-stage local fdr procedures and Storey's procedure.
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R")
#          eb_rejections(x, level = 0.05)   # number of positions declared non-null
# Source:  simEM_fcns.R (gpd, gpdMixture, tau0 to tau3, f0 to f14, lfdr, Dn,
#          zigpEMAlgorithm, ModifiedBH, ModifiedBHFDR, negLogLh), simFDR_fcns.R
#          (EM_fdr) and real_scenarios.R (RealDataFDR).
# -----------------------------------------------------------------------------

# Small constant the fit adds to probabilities and to denominators.
EB_FLOOR <- 10^(-10)

# ----- ZIGP null fitted by EM -----

#' Generalized Poisson probabilities.
#' @param x Counts.
#' @param lambda,theta Rate and dispersion.
#' @return Vector of probabilities.
gp_pmf <- function(x, lambda, theta) {
  exp(log(lambda) + (x - 1) * log(lambda + theta * x) - lfactorial(x) - lambda - theta * x)
}

#' ZIGP probabilities, kept away from zero.
#' @param x Counts.
#' @param eta,lambda,theta ZIGP parameters.
#' @return Vector of probabilities.
zigp_pmf_floored <- function(x, eta, lambda, theta) {
  eta * (1 * (x == 0)) + (1 - eta) * (gp_pmf(x, lambda, theta) + EB_FLOOR) + EB_FLOOR
}

#' Probability that a ZIGP count is at most C.
#' @param eta,lambda,theta ZIGP parameters.
#' @param C Cutoff.
#' @return The probability.
zigp_cdf <- function(eta, lambda, theta, C) {
  eta + (1 - eta) * sum(gp_pmf(0:C, lambda, theta))
}

#' Expected contributions of the counts above the cutoff in the E-step.
#'
#' The fit sees only the n counts up to C; the counts above C are replaced by
#' their expectation under the current ZIGP, summed up to the largest count K.
#'
#' @param K Largest observed count.
#' @param eta,lambda,theta Current ZIGP parameters.
#' @param C Cutoff.
#' @param n Number of counts up to C.
#' @return List of the expected number of counts above C (`n`), their sum
#'   (`total`), and the sums weighted by lambda / (lambda + theta x) and by
#'   theta x / (lambda + theta x) (`total_rate`, `n_rate`, `total_dispersion`,
#'   `n_dispersion`).
em_tail <- function(K, eta, lambda, theta, C, n) {
  below <- zigp_cdf(eta, lambda, theta, C)
  v <- (C + 1):K
  pmf <- gp_pmf(v, lambda, theta)
  pmf_times_v <- exp(log(v * lambda) + (v - 1) * log(lambda + theta * v) - lfactorial(v) -
                       lambda - theta * v)
  rate <- lambda / (lambda + theta * v)
  dispersion <- (theta * v) / (lambda + theta * v)
  scaled <- function(terms) n * ((1 - eta) * sum(terms)) / below
  list(n = n * ((1 - eta) * (1 - sum(gp_pmf(0:C, lambda, theta)))) / below,
       total = scaled(pmf_times_v),
       total_rate = scaled(pmf_times_v * rate), n_rate = scaled(pmf * rate),
       total_dispersion = scaled(pmf_times_v * dispersion), n_dispersion = scaled(pmf * dispersion))
}

#' EM estimates of the ZIGP null for every cutoff C = 1, ..., max(x).
#'
#' For each C the ZIGP is fitted to the counts up to C, treating the counts
#' above C as missing. Starting values are moment estimates from the lower
#' three quarters of the data; the iterations stop when no parameter moves
#' by more than 1e-4.
#'
#' @param x Counts.
#' @return Matrix with rows eta, lambda, theta and one column per cutoff.
zigp_em <- function(x) {
  K <- max(x)
  epsilon <- 10^(-4)
  estimates <- matrix(0, nrow = 3, ncol = K, dimnames = list(c("eta", "lambda", "theta"), NULL))

  lower <- x[which(x <= quantile(x, 0.75))]
  ex <- mean(lower)
  vx <- var(lower)
  er <- mean(x == 0) / (1 - mean(x == 0))
  phi <- sqrt(vx / ex - er * ex + 10^(-10))
  if (ex == 0 || phi == 0 || is.nan(phi)) {
    theta_start <- 1
    lambda_start <- 0
  } else {
    theta_start <- min(max(0, 1 - 1 / phi), 1)
    lambda_start <- max(0, ex / ((1 - mean(x == 0)) * phi))
  }
  eta_start <- max(0, mean(x == 0))

  for (C in 1:K) {
    count <- c(length(x[which(x == 0)]), tabulate(x, nbins = C))
    n <- sum(count)
    values <- 0:C
    lagged <- c(0, seq(0, C - 1, 1))
    eta <- eta_start
    lambda <- lambda_start
    theta <- theta_start
    repeat {
      eta0 <- max(0, eta)
      lambda0 <- max(EB_FLOOR, lambda)
      theta0 <- max(EB_FLOOR, theta)
      tail <- em_tail(K, eta0, lambda0, theta0, C, n)
      mixture <- zigp_pmf_floored(values, eta0, lambda0, theta0)
      from_gp <- (1 - eta0) * (gp_pmf(values, lambda0, theta0) + EB_FLOOR) / mixture
      u1 <- count[1] * (eta0 * (1 * (0 == 0)) / zigp_pmf_floored(0, eta0, lambda0, theta0))
      u2 <- sum(count * from_gp)
      eta <- u1 / (u1 + u2 + tail$n + EB_FLOOR)
      A <- sum(values * count * from_gp) + tail$total
      u6 <- sum(lagged * count * (lambda0 / (lambda0 + theta0 * values)))
      B <- u2 + tail$n + u6 + tail$total_rate - tail$n_rate
      lambda <- (B + EB_FLOOR) / (u2 + tail$n + EB_FLOOR)
      v1 <- sum(lagged * count * ((theta0 * values) / (lambda0 + theta0 * values)))
      D <- v1 + tail$total_dispersion - tail$n_dispersion
      theta <- (D + EB_FLOOR) / (A + EB_FLOOR)
      moved <- abs(eta0 - eta) > epsilon || abs(lambda0 - lambda) > epsilon ||
        abs(theta0 - theta) > epsilon
      if (!moved) break
    }
    estimates[, C] <- c(eta, lambda, theta)
  }
  estimates
}

# ----- Cutoff and null proportion -----

#' Negative log likelihood of the ZIGP over the counts up to C.
#' @param count Number of observations at 0, 1, ..., max(x).
#' @param eta,lambda,theta ZIGP parameters.
#' @param C Cutoff.
#' @return The negative log likelihood.
zigp_negloglik <- function(count, eta, lambda, theta, C) {
  v <- 1:C
  at_zero <- count[1] * (-log(eta + (1 - eta) * exp(-lambda)))
  positive <- count[2:(C + 1)] * (-log((1 - eta) * exp(log(lambda) + (v - 1) * log(lambda + theta * v) -
                                                         lfactorial(v) - lambda - theta * v)))
  at_zero + sum(positive)
}

#' Fit the ZIGP null, the null proportion and the cutoff by empirical Bayes.
#'
#' For each candidate cutoff C the ZIGP is fitted by EM and pi0 is the share
#' of counts up to C divided by the fitted probability of that range. The
#' cutoff is the first local minimum (or else the minimum) of the negative
#' log likelihood of the counts up to C plus that of the counts above C under
#' their observed frequencies.
#'
#' @param x Counts.
#' @return List: `estimates` (matrix from `zigp_em()`), `pi0` (per cutoff),
#'   `criterion` (per cutoff), `cutoff`, and `eta`, `lambda`, `theta`, `pi0_hat`
#'   at the chosen cutoff.
eb_fit <- function(x) {
  K <- max(x)
  M <- length(x)
  count <- c(sum(x == 0), tabulate(x, K))
  estimates <- zigp_em(x)
  pi0 <- numeric(K)
  criterion <- numeric(K)
  for (C in 1:K) {
    eta <- estimates[1, C]
    lambda <- estimates[2, C]
    theta <- estimates[3, C]
    pi0[C] <- min(1, sum(count[1:(C + 1)]) / (M * zigp_cdf(eta, lambda, theta, C)))
    above <- if (C < K) {
      log(M) * length(x[which(x > C)]) -
        sum(count[(C + 2):(K + 1)] * (log(count[(C + 2):(K + 1)] + EB_FLOOR)))
    } else {
      log(M) * length(x[which(x > C)]) - sum(count[K] * (log(count[K] + EB_FLOOR)))
    }
    criterion[C] <- above + zigp_negloglik(count, eta, lambda, theta, C)
  }
  cutoff <- min(which(diff(sign(diff(criterion))) > 0) + 2, which.min(criterion))
  list(estimates = estimates, pi0 = pi0, criterion = criterion, cutoff = cutoff,
       eta = estimates[1, cutoff], lambda = estimates[2, cutoff], theta = estimates[3, cutoff],
       pi0_hat = pi0[cutoff])
}

#' Empirical Bayes local fdr of sorted counts.
#'
#' fdr(x) = pi0 f0(x) / f(x), with f(x) the observed relative frequency of x.
#'
#' @param fit Output of `eb_fit()`.
#' @param x Counts at which to evaluate.
#' @param inverse_frequency Number of observations divided by the number of
#'   them equal to each element of `x`.
#' @return Vector of local fdr values.
eb_local_fdr <- function(fit, x, inverse_frequency) {
  fit$pi0_hat * (fit$eta * (1 * (x == 0)) + (1 - fit$eta) * (gp_pmf(x, fit$lambda, fit$theta) + EB_FLOOR) +
                   EB_FLOOR) * inverse_frequency
}

#' Number of positions declared non-null in a data set (one-stage procedure).
#'
#' @param x Counts, one per position.
#' @param level fdr level.
#' @return Number of positions with local fdr below `level`.
eb_rejections <- function(x, level = 0.05) {
  fit <- eb_fit(x)
  sorted <- sort(x)
  values <- sort(unique(x))
  frequency <- as.numeric(table(factor(x, levels = values)))
  fdr <- eb_local_fdr(fit, sorted, length(x) / frequency[match(sorted, values)])
  sum(fdr[!is.na(fdr)] < level)
}

# ----- The three procedures in a simulated data set -----

#' Second-stage threshold: counts from here on are rejected outright.
#' @param fit Output of `eb_fit()`.
#' @param n Number of counts.
#' @return The threshold.
second_stage_threshold <- function(fit, n) {
  bound <- max(fit$lambda / (exp(fit$theta - 1) - fit$theta), log(n) / (fit$theta - 1 - log(fit$theta)))
  ceiling(max(fit$cutoff + 1, bound))
}

#' Rejections of Storey's procedure with ZIGP tail p-values.
#'
#' The p-value of a count is its ZIGP upper tail probability. The sorted
#' p-values are compared with (number of counts at least as large as each
#' distinct value) * level / (n pi0).
#'
#' @param data Output of `simulated_data()`.
#' @param fit Output of `eb_fit()`.
#' @param level FDR level.
#' @return Positions (among the distinct counts) that are rejected.
storey_rejected <- function(data, fit, level) {
  n <- length(data$x)
  largest <- max(data$x)
  tail <- vapply(seq_len(largest), function(i) {
    max(0, (1 - fit$eta) * (1 - sum(gp_pmf(0:(i - 1), fit$lambda, fit$theta))))
  }, numeric(1))
  pvalue <- c(1, tail)[(0:largest) %in% data$values]
  frequency <- data$null_frequency + data$alt_frequency

  U <- length(pvalue)
  by_pvalue <- order(pvalue)
  at_least <- numeric()
  for (i in 2:length(frequency)) {
    at_least[1] <- n
    at_least[i] <- n - sum(frequency[1:(i - 1)])
  }
  below <- pvalue[by_pvalue] <= at_least * (level / (n * fit$pi0_hat))
  if (sum(below) == 0) return(integer(0))
  last <- max(which(below))
  if (last == 1) return(integer(0))  # a single rejection is not kept
  by_pvalue[1:last]
}

#' Rejection counts of the two-stage, one-stage and Storey procedures in one
#' simulated data set.
#'
#' Two-stage: null counts below the second-stage threshold are rejected by
#' their local fdr, every count from the threshold on is rejected. One-stage:
#' counts are rejected by their local fdr alone.
#'
#' @param data Output of `simulated_data()`.
#' @param C_true Cutoff the data were simulated with.
#' @param level fdr level.
#' @return One-row data frame with TS_R, TS_V, TS_TP, OS_R, OS_V, OS_TP,
#'   Storey_R, Storey_V, Storey_TP (rejections, false and true rejections).
eb_simulation_counts <- function(data, C_true, level) {
  null <- data$null
  alt <- data$alt
  n <- length(data$x)
  fit <- eb_fit(data$x)
  threshold <- second_stage_threshold(fit, n)

  null_sorted <- sort(null)
  alt_sorted <- sort(alt)
  weight <- function(x) n / (data$null_frequency + data$alt_frequency)[match(x, data$values)]
  alt_weight <- weight(alt_sorted)
  below_level <- function(x, weights = weight(x)) {
    fdr <- eb_local_fdr(fit, x, weights)
    sum(fdr[!is.na(fdr)] < level)
  }

  # Null counts rejected by their local fdr.
  null_below_threshold <- below_level(null_sorted[null_sorted < threshold])
  null_up_to_cutoff <- below_level(null_sorted[null_sorted <= fit$cutoff])
  null_above_cutoff <- below_level(null_sorted[null_sorted > fit$cutoff])
  # Non-null counts rejected by their local fdr; the weights are taken from
  # the smallest non-null counts in both ranges.
  n_alt_below <- sum(alt < threshold)
  n_alt_from <- sum(alt > threshold - 1)
  alt_below_threshold <- below_level(alt_sorted[alt_sorted < threshold], alt_weight[seq_len(n_alt_below)])
  alt_from_threshold <- below_level(alt_sorted[alt_sorted > (threshold - 1)], alt_weight[seq_len(n_alt_from)])

  TS_V <- max(0, null_below_threshold - null_up_to_cutoff) + sum(null >= threshold)
  TS_TP <- alt_below_threshold + sum(alt >= threshold)
  OS_V <- max(0, null_above_cutoff)
  OS_TP <- alt_below_threshold + alt_from_threshold

  storey <- storey_rejected(data, fit, level)
  Storey_V <- sum(data$null_frequency[storey])
  Storey_TP <- sum(data$alt_frequency[storey])

  data.frame(TS_R = TS_TP + TS_V, TS_V = TS_V, TS_TP = TS_TP,
             OS_R = OS_TP + OS_V, OS_V = OS_V, OS_TP = OS_TP,
             Storey_R = Storey_V + Storey_TP, Storey_V = Storey_V, Storey_TP = Storey_TP)
}
