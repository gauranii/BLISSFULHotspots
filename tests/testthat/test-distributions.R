# -----------------------------------------------------------------------------
# File:    tests/testthat/test-distributions.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Unit tests of R/distributions.R: the ZIGP probabilities against
#          the Poisson and their sum, the three samplers against their means
#          and under a seed, and the shifted alternative densities.
# Inputs:  None.
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/setup.R")

test_that("dzigp sums to one and has the ZIGP mean", {
  x <- 0:400
  p <- dzigp(x, eta = 0.3, lambda = 2, theta = 0.2)
  expect_length(p, length(x))
  expect_true(all(p >= 0))
  expect_equal(sum(p), 1, tolerance = 1e-10)
  expect_equal(sum(x * p), (1 - 0.3) * 2 / (1 - 0.2), tolerance = 1e-10)
})

test_that("dzigp with theta = 0 is the zero-inflated Poisson", {
  x <- c(3, 0, 7, 0, 1)
  expect_equal(dzigp(x, 0.4, 1.5, 0), 0.4 * (x == 0) + 0.6 * dpois(x, 1.5))
  expect_equal(dzigp(0:5, 0, 2.5, 0), dpois(0:5, 2.5))
  # worked by hand: P(0) = eta + (1 - eta) exp(-lambda), P(1) = (1 - eta) lambda exp(-lambda - theta)
  expect_equal(dzigp(c(0, 1), 0.25, 2, 0.5), c(0.25 + 0.75 * exp(-2), 0.75 * 2 * exp(-2.5)))
})

test_that("dzigp checks the range of the parameters unless told not to", {
  expect_error(dzigp(1, 1.2, 2, 0.1))
  expect_error(dzigp(1, 0.2, 0, 0.1))
  expect_error(dzigp(1, 0.2, 2, 1.1))
  expect_equal(dzigp(c(0, 2), 0.2, 2, 1.1, check = FALSE),
               c(0.2 + 0.8 * exp(-2), 0.8 * 2 * (2 + 2.2) * exp(-2 - 2.2) / 2))
})

test_that("rzigp gives the same draw under the same seed and uses one uniform number", {
  set.seed(31)
  first <- vapply(1:20, function(i) rzigp(0.4, 2, 0.3), numeric(1))
  after_first <- runif(1)
  set.seed(31)
  second <- vapply(1:20, function(i) rzigp(0.4, 2, 0.3), numeric(1))
  expect_equal(first, second)
  set.seed(31)
  uniforms <- runif(21)
  expect_equal(after_first, uniforms[21])
})

test_that("rzigp inverts the distribution function", {
  # with eta = 0 and theta = 0 the draw is the Poisson quantile of the uniform number
  set.seed(7)
  draws <- vapply(1:300, function(i) rzigp(0, 3, 0), numeric(1))
  set.seed(7)
  expect_equal(draws, qpois(runif(300), 3))
  # a uniform number below P(X = 0) gives zero
  expect_equal(vapply(1:25, function(i) rzigp(1, 2, 0.2), numeric(1)), rep(0, 25))
})

test_that("rzigp has the ZIGP mean and share of zeros over many draws", {
  set.seed(102)
  draws <- vapply(1:6000, function(i) rzigp(0.4, 2, 0.3), numeric(1))
  expect_equal(mean(draws), 0.6 * 2 / 0.7, tolerance = 0.05)
  expect_equal(mean(draws == 0), 0.4 + 0.6 * exp(-2), tolerance = 0.05)
  expect_equal(as.numeric(table(factor(draws, 0:3))) / 6000, dzigp(0:3, 0.4, 2, 0.3),
               tolerance = 0.1)
  expect_error(rzigp(0.4, -1, 0.3))
})

test_that("rdirichlet is the normalized gamma draws and sums to one", {
  beta <- c(0.5, 2, 1, 4)
  set.seed(5)
  draw <- rdirichlet(beta)
  set.seed(5)
  gammas <- rgamma(4, beta, 1)
  expect_equal(draw, gammas / sum(gammas))
  expect_equal(sum(draw), 1)
  set.seed(5)
  expect_identical(rdirichlet(beta), draw)
})

test_that("rdirichlet has mean beta / sum(beta) over many draws", {
  beta <- c(1, 2, 5)
  set.seed(77)
  draws <- t(replicate(4000, rdirichlet(beta)))
  expect_equal(colMeans(draws), beta / sum(beta), tolerance = 0.03)
})

test_that("rdirichlet replaces a missing concentration only when asked", {
  set.seed(1)
  expect_true(suppressWarnings(anyNA(rdirichlet(c(1, NA, 2)))))
  set.seed(1)
  draw <- rdirichlet(c(1, NA, 2), replace_missing = TRUE)
  expect_false(anyNA(draw))
  expect_equal(sum(draw), 1)
  expect_lt(draw[2], 0.01)  # concentration 1e-5
})

test_that("rmultinomial keeps the number of trials and skips empty categories", {
  set.seed(9)
  draw <- rmultinomial(50, c(0.2, 0, 0.5, 0.3))
  expect_equal(sum(draw), 50)
  expect_equal(draw[2], 0)
  expect_equal(rmultinomial(7, c(0, 1, 0)), c(0, 7, 0))
  set.seed(9)
  expect_identical(rmultinomial(50, c(0.2, 0, 0.5, 0.3)), draw)
})

test_that("rmultinomial has mean size * prob, also for weights that do not sum to one", {
  set.seed(12)
  draws <- t(replicate(3000, rmultinomial(20, c(1, 3, 4, 2))))
  expect_true(all(rowSums(draws) == 20))
  expect_equal(colMeans(draws), 20 * c(0.1, 0.3, 0.4, 0.2), tolerance = 0.03)
  # one trial: an indicator of the category drawn
  one <- t(replicate(2000, rmultinomial(1, c(0.25, 0.75))))
  expect_true(all(rowSums(one) == 1))
  expect_equal(mean(one[, 2]), 0.75, tolerance = 0.05)
})

test_that("shifted_density is a Poisson that starts right above the cutoff", {
  x <- 0:60
  f1 <- shifted_density("P", x, xi = 3, delta = NA, nu = NA, C = 5)
  expect_equal(f1[x <= 5], rep(0, 6))
  expect_equal(f1[x > 5], dpois(x[x > 5] - 6, 3))
  expect_equal(sum(f1), 1, tolerance = 1e-10)
})

test_that("shifted_density is a generalized Poisson that starts right above the cutoff", {
  x <- 0:300
  f1 <- shifted_density("GP", x, xi = NA, delta = 2, nu = 0.3, C = 4)
  expect_equal(f1[x <= 4], rep(0, 5))
  expect_equal(f1[x > 4], dzigp(x[x > 4] - 5, 0, 2, 0.3))
  expect_equal(sum(f1), 1, tolerance = 1e-10)
  # nu = 0: a Poisson with mean delta
  expect_equal(shifted_density("GP", 0:12, NA, 2.5, 0, 2)[4:13], dpois(0:9, 2.5))
  # another model has no parametric f1
  expect_equal(shifted_density("NP", 0:5, 1, 1, 0.5, 2), rep(0, 6))
})
