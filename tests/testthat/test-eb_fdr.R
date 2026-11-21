# -----------------------------------------------------------------------------
# File:    tests/testthat/test-eb_fdr.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Unit tests of R/eb_fdr.R: the ZIGP probabilities and the E-step
#          sums against their definitions, the EM fit on counts drawn from
#          known parameters, the cutoff and null proportion on a simulated
#          mixture, and the local fdr, the thresholds and the rejections of
#          the three procedures on cases worked by hand.
# Inputs:  None.
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/setup.R")

# A mixture with the two groups far apart: null counts 0 to 6, non-null 11 up.
SEPARATED <- simulated_data(small_setting(N = 400), seed = 11)
SEPARATED_FIT <- eb_fit(SEPARATED$x)
# A fit with round numbers, for the cases worked by hand.
HAND_FIT <- list(eta = 0.4, lambda = 1.5, theta = 0, pi0_hat = 0.9, cutoff = 5)

test_that("gp_pmf sums to one and is the Poisson when theta = 0", {
  expect_equal(gp_pmf(0:8, 2.5, 0), dpois(0:8, 2.5))
  expect_equal(sum(gp_pmf(0:400, 2, 0.3)), 1, tolerance = 1e-10)
  expect_equal(sum(0:400 * gp_pmf(0:400, 2, 0.3)), 2 / 0.7, tolerance = 1e-10)
  expect_equal(gp_pmf(0:6, 2, 0.3), dzigp(0:6, 0, 2, 0.3))
})

test_that("zigp_pmf_floored is the ZIGP probability plus the two floors", {
  x <- 0:10
  expect_equal(zigp_pmf_floored(x, 0.4, 2, 0.3) - dzigp(x, 0.4, 2, 0.3),
               rep(0.6 * EB_FLOOR + EB_FLOOR, 11), tolerance = 1e-6)
  expect_true(all(zigp_pmf_floored(c(200, 500), 0.4, 2, 0.3) > 0))
})

test_that("zigp_cdf sums the ZIGP probabilities up to the cutoff", {
  expect_equal(zigp_cdf(0.4, 2, 0.3, 5), sum(dzigp(0:5, 0.4, 2, 0.3)))
  expect_equal(zigp_cdf(0.4, 2, 0.3, 0), 0.4 + 0.6 * exp(-2))
  expect_equal(zigp_cdf(0.4, 2, 0.3, 400), 1, tolerance = 1e-10)
})

test_that("em_tail gives the expected counts above the cutoff and their sums", {
  eta <- 0.4; lambda <- 2; theta <- 0.3; C <- 4; n <- 90
  tail <- em_tail(K = 400, eta, lambda, theta, C, n)
  below <- sum(dzigp(0:C, eta, lambda, theta))
  v <- 5:400
  expect_named(tail, c("n", "total", "total_rate", "n_rate", "total_dispersion", "n_dispersion"))
  # n counts up to C stand for n (1 - F(C)) / F(C) counts above it
  expect_equal(tail$n, n * (1 - below) / below)
  expect_equal(tail$total, n * sum(v * dzigp(v, eta, lambda, theta)) / below)
  # the two weights lambda / (lambda + theta x) and theta x / (lambda + theta x) sum to one
  expect_equal(tail$total_rate + tail$total_dispersion, tail$total)
  expect_equal(tail$n_rate + tail$n_dispersion, tail$n, tolerance = 1e-8)
  expect_equal(tail$n_rate, n * sum(dzigp(v, eta, lambda, theta) * lambda / (lambda + theta * v)) / below)
})

test_that("zigp_em recovers the parameters of zero-inflated Poisson counts", {
  set.seed(3)
  x <- vapply(1:4000, function(i) rzigp(0.4, 2, 0), numeric(1))
  estimates <- zigp_em(x)
  expect_equal(dim(estimates), c(3L, max(x)))
  expect_equal(rownames(estimates), c("eta", "lambda", "theta"))
  at_largest <- estimates[, max(x)]
  expect_equal(at_largest[["eta"]], 0.4, tolerance = 0.1)
  expect_equal(at_largest[["lambda"]], 2, tolerance = 0.1)
  expect_lt(at_largest[["theta"]], 0.05)
  expect_true(all(estimates["eta", ] >= 0 & estimates["eta", ] <= 1))
  expect_true(all(estimates["lambda", ] > 0))
  expect_identical(zigp_em(x), estimates)  # no random number is used
})

test_that("zigp_negloglik is minus the ZIGP log likelihood of the counts up to C", {
  count <- c(40, 25, 12, 6, 3, 1)
  expect_equal(zigp_negloglik(count, 0.3, 1.8, 0.2, C = 5),
               -sum(count * log(dzigp(0:5, 0.3, 1.8, 0.2))))
  expect_equal(zigp_negloglik(count, 0.3, 1.8, 0.2, C = 2),
               -sum(count[1:3] * log(dzigp(0:2, 0.3, 1.8, 0.2))))
})

test_that("eb_fit finds the cutoff and the null proportion of a separated mixture", {
  fit <- SEPARATED_FIT
  K <- max(SEPARATED$x)
  expect_named(fit, c("estimates", "pi0", "criterion", "cutoff", "eta", "lambda", "theta", "pi0_hat"))
  expect_equal(dim(fit$estimates), c(3L, K))
  expect_length(fit$pi0, K)
  expect_length(fit$criterion, K)
  expect_true(all(fit$pi0 > 0 & fit$pi0 <= 1))
  # the data have null counts up to 6 and non-null counts from 11: C = 10 separates them
  expect_gte(fit$cutoff, max(SEPARATED$null))
  expect_lt(fit$cutoff, min(SEPARATED$alt))
  expect_equal(fit$pi0_hat, length(SEPARATED$null) / 400, tolerance = 0.05)
  expect_equal(fit$eta[["eta"]], 0.4, tolerance = 0.2)
  expect_equal(c(fit$eta, fit$lambda, fit$theta), fit$estimates[, fit$cutoff])
  expect_equal(fit$pi0_hat, fit$pi0[fit$cutoff])
  # here the criterion is smallest at the cutoff
  expect_equal(fit$criterion[fit$cutoff], min(fit$criterion))
})

test_that("eb_fit takes the minimum of a criterion without a local minimum", {
  x <- c(rep(0, 12), rep(1, 6), 2, 2)
  fit <- eb_fit(x)
  expect_length(fit$criterion, 2)
  expect_equal(fit$cutoff, which.min(fit$criterion))
})

test_that("eb_local_fdr is pi0 f0(x) over the relative frequency of x", {
  x <- c(0, 2, 7)
  inverse_frequency <- c(2, 10, 50)
  expect_equal(eb_local_fdr(HAND_FIT, x, inverse_frequency),
               0.9 * (dzigp(x, 0.4, 1.5, 0) + 0.6 * EB_FLOOR + EB_FLOOR) * inverse_frequency)
  # a count that all of 20 observations share, under a point mass at zero
  point_mass <- list(eta = 1, lambda = 1, theta = 0, pi0_hat = 1)
  expect_equal(eb_local_fdr(point_mass, 0, 1), 1, tolerance = 1e-8)
})

test_that("eb_rejections counts the positions with local fdr below the level", {
  x <- SEPARATED$x
  expect_equal(eb_rejections(x, 0.05), length(SEPARATED$alt))
  # the same count from the fit, one distinct count at a time
  values <- sort(unique(x))
  frequency <- as.numeric(table(x))
  fdr <- eb_local_fdr(SEPARATED_FIT, values, length(x) / frequency)
  expect_equal(eb_rejections(x, 0.05), sum(frequency[fdr < 0.05]))
  expect_equal(eb_rejections(x, 0), 0)
  expect_equal(eb_rejections(sample(x)), eb_rejections(x))  # the order of the positions is free
})

test_that("second_stage_threshold is worked by hand", {
  # lambda / (exp(theta - 1) - theta) = 2 / (exp(-0.8) - 0.2) = 8.02 and
  # log(n) / (theta - 1 - log(theta)) = log(100) / (log(5) - 0.8) = 5.69
  fit <- list(lambda = 2, theta = 0.2, cutoff = 5)
  expect_equal(2 / (exp(-0.8) - 0.2), 8.02, tolerance = 1e-3)
  expect_equal(second_stage_threshold(fit, 100), 9)
  # never below the count right above the cutoff
  expect_equal(second_stage_threshold(modifyList(fit, list(cutoff = 20)), 100), 21)
  # more counts raise the second bound: log(1e6) / (log(5) - 0.8) = 17.07
  expect_equal(second_stage_threshold(fit, 1e6), 18)
})

test_that("storey_rejected is worked by hand", {
  # p-values of the counts 0, 1, 2 under the fit: 1, 0.6 (1 - P(0)) = 0.466 and
  # 0.6 (1 - P(0) - P(1)) = 0.265; thresholds (100, 50, 20) * 0.05 / (100 * 0.9): none below
  data <- list(x = rep(c(0, 1, 2), c(50, 30, 20)), values = c(0, 1, 2),
               null_frequency = c(50, 30, 20), alt_frequency = c(0, 0, 0))
  expect_equal(0.6 * (1 - dpois(0, 1.5)), 0.466, tolerance = 1e-3)
  expect_equal(storey_rejected(data, HAND_FIT, 0.05), integer(0))
  # two large counts with p-values near zero are rejected
  data <- list(x = rep(c(0, 1, 30, 31), c(60, 30, 5, 5)), values = c(0, 1, 30, 31),
               null_frequency = c(60, 30, 0, 0), alt_frequency = c(0, 0, 5, 5))
  expect_equal(sort(storey_rejected(data, HAND_FIT, 0.05)), 3:4)
  # a single rejection is not kept
  data <- list(x = rep(c(0, 1, 30), c(60, 30, 10)), values = c(0, 1, 30),
               null_frequency = c(60, 30, 0), alt_frequency = c(0, 0, 10))
  expect_equal(storey_rejected(data, HAND_FIT, 0.05), integer(0))
})

test_that("storey_rejected rejects the non-null counts of a separated mixture", {
  rejected <- storey_rejected(SEPARATED, SEPARATED_FIT, 0.05)
  expect_equal(sum(SEPARATED$alt_frequency[rejected]), length(SEPARATED$alt))
  expect_true(all(SEPARATED$values[rejected] > max(SEPARATED$null)))
})

test_that("eb_simulation_counts finds every non-null count of a separated mixture", {
  counts <- eb_simulation_counts(SEPARATED, C_true = 10, level = 0.05)
  n_alt <- length(SEPARATED$alt)
  expect_equal(dim(counts), c(1L, 9L))
  expect_named(counts, c("TS_R", "TS_V", "TS_TP", "OS_R", "OS_V", "OS_TP",
                         "Storey_R", "Storey_V", "Storey_TP"))
  expect_equal(unlist(counts[c("TS_TP", "OS_TP", "Storey_TP")], use.names = FALSE), rep(n_alt, 3))
  expect_equal(unlist(counts[c("TS_V", "OS_V", "Storey_V")], use.names = FALSE), rep(0, 3))
  expect_equal(counts$OS_R, eb_rejections(SEPARATED$x, 0.05))
})

test_that("eb_simulation_counts adds up: R = V + TP within the group sizes", {
  data <- simulated_data(small_setting(N = 500, lambda = 2, theta = 0.3, C = 5), seed = 5)
  counts <- eb_simulation_counts(data, C_true = 5, level = 0.05)
  for (procedure in c("TS", "OS", "Storey")) {
    R <- counts[[paste0(procedure, "_R")]]
    V <- counts[[paste0(procedure, "_V")]]
    TP <- counts[[paste0(procedure, "_TP")]]
    expect_equal(R, V + TP, info = procedure)
    expect_true(V >= 0 && V <= length(data$null), info = procedure)
    expect_true(TP >= 0 && TP <= length(data$alt), info = procedure)
  }
  expect_gt(counts$TS_TP, 0.5 * length(data$alt))
  expect_identical(eb_simulation_counts(data, 5, 0.05), counts)
})
