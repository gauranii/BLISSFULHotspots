# -----------------------------------------------------------------------------
# File:    tests/testthat/test-bayes_fdr.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Unit tests of R/bayes_fdr.R: the mixture density and local fdr on
#          a case worked by hand, each full conditional draw against its
#          distribution and under a seed, the log likelihoods and log
#          posteriors against their definitions, the proposal covariance,
#          and short runs of the sampler for the four alternatives.
# Inputs:  None.
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/setup.R")

V2018 <- SAMPLER_VERSIONS[["2018-03-25"]]
V2021 <- SAMPLER_VERSIONS[["2021-02-19"]]

# Distinct counts 0 to 8 with their frequencies, and a state of the chain.
U_VALUES <- 0:8
FREQUENCY <- c(20, 10, 6, 3, 1, 2, 3, 3, 2)
X_SORTED <- rep(U_VALUES, FREQUENCY)
STATE <- list(eta = 0.4, lambda = 1.5, theta = 0.1, pi0 = 0.8, xi = 2, delta = 2, nu = 0.2, tau = 3,
              psi_sp = rep(1 / 9, 9), psi_np = rep(1 / 9, 9))

# A data set with the two groups far apart, for the runs of the sampler.
SEPARATED <- simulated_data(small_setting(N = 60), seed = 21)

test_that("SAMPLER_VERSIONS lists the three versions and what differs between them", {
  expect_equal(BAYES_MODELS, c("P", "GP", "SP", "NP"))
  expect_named(SAMPLER_VERSIONS, c("2018-03-25", "2021-02-07", "2021-02-19"))
  for (version in SAMPLER_VERSIONS) {
    expect_named(version, c("check_range", "replace_missing", "no_cutoff", "extra_atoms",
                            "logit_scale", "theta_max", "sp_prior", "fdr_from"))
  }
  expect_equal(V2018$extra_atoms(U = 60, largest = 80), 120)   # 1.5 times the largest count
  expect_equal(V2021$extra_atoms(U = 60, largest = 80), 70)    # 10 more than the distinct counts
  expect_equal(V2018$fdr_from(burn = 1000), 1)
  expect_equal(V2021$fdr_from(burn = 1000), 1001)
  expect_equal(c(V2018$theta_max, V2021$theta_max), c(1, 0.5))
})

test_that("sampler_settings holds the starting values and the priors", {
  settings <- sampler_settings()
  expect_named(settings, c("init", "pi0_prior", "tau_prior", "xi_prior", "proposal", "adapt"))
  expect_equal(c(settings$init$lambda, settings$init$theta), c(2, 0.25))
  expect_equal(settings$init$atoms, 50)
  expect_equal(sum(settings$init$beta), 1)
  expect_equal(settings$pi0_prior, list(alpha = 1, beta = 1))
  expect_equal(settings$proposal$null_var, diag(0.1, 3))
  expect_equal(settings$proposal$gp_var, diag(0.1, 2))
  other <- sampler_settings(lambda_init = 1, theta_init = 0.5)
  expect_equal(c(other$init$lambda, other$init$theta), c(1, 0.5))
})

test_that("mixture_fdr is pi0 f0 / f, and one up to the cutoff (Poisson f1)", {
  out <- mixture_fdr("P", U_VALUES, STATE, C = 3, V2021)
  f0 <- 0.8 * dzigp(U_VALUES, 0.4, 1.5, 0.1)
  f1 <- c(0, 0, 0, 0, 0.2 * dpois(0:4, 2))
  expect_equal(out$f0, f0)
  expect_equal(out$f, f0 + f1)
  expect_equal(out$fdr, f0 / (f0 + f1))
  expect_equal(out$fdr[1:4], rep(1, 4))
  expect_true(all(out$fdr[5:9] < 1))
})

test_that("mixture_fdr uses the f1 of each model above the cutoff", {
  f0 <- 0.8 * dzigp(U_VALUES, 0.4, 1.5, 0.1)
  above <- U_VALUES > 3
  gp <- mixture_fdr("GP", U_VALUES, STATE, 3, V2021)
  expect_equal(gp$f[above], f0[above] + 0.2 * dzigp(0:4, 0, 2, 0.2))
  sp <- mixture_fdr("SP", U_VALUES, STATE, 3, V2021)
  expect_equal(sp$f[above], f0[above] + 0.2 / 9)
  np <- mixture_fdr("NP", U_VALUES, STATE, 3, V2021)
  expect_equal(np$f[above], rep(1 / 9, 5))
  expect_equal(np$fdr[above], pmin(1, f0[above] / (1 / 9)))
  for (out in list(gp, sp, np)) {
    expect_equal(out$f[!above], f0[!above])
    expect_true(all(out$fdr >= 0 & out$fdr <= 1))
  }
  # f0 given: used as it is
  expect_equal(mixture_fdr("NP", U_VALUES, STATE, 3, V2021, f0 = rep(0.05, 9))$fdr[above], rep(0.45, 5))
})

test_that("mixture_fdr without a cutoff uses C = 1 or stops, by version", {
  expect_equal(mixture_fdr("P", U_VALUES, STATE, NULL, V2021), mixture_fdr("P", U_VALUES, STATE, 1, V2021))
  expect_error(mixture_fdr("P", U_VALUES, STATE, NULL, V2018), "no cutoff could be drawn")
})

test_that("draw_cutoff draws one of the distinct counts, the same under the same seed", {
  tau_prior <- sampler_settings()$tau_prior
  set.seed(4)
  draws <- vapply(1:200, function(i) draw_cutoff("P", U_VALUES, FREQUENCY, STATE, tau_prior, V2021),
                  numeric(1))
  expect_true(all(draws %in% U_VALUES))
  set.seed(4)
  again <- vapply(1:200, function(i) draw_cutoff("P", U_VALUES, FREQUENCY, STATE, tau_prior, V2021),
                  numeric(1))
  expect_equal(draws, again)
})

test_that("draw_cutoff follows its full conditional", {
  tau_prior <- sampler_settings()$tau_prior
  # log posterior of each candidate: the log likelihood of the mixture with that cutoff plus the
  # Poisson(tau) log prior
  logpost <- vapply(U_VALUES, function(C) {
    sum(FREQUENCY * log(mixture_fdr("P", U_VALUES, STATE, C, V2021)$f)) + dpois(C, STATE$tau, log = TRUE)
  }, numeric(1))
  probability <- exp(logpost - max(logpost)) / sum(exp(logpost - max(logpost)))
  set.seed(8)
  draws <- vapply(1:3000, function(i) draw_cutoff("P", U_VALUES, FREQUENCY, STATE, tau_prior, V2021),
                  numeric(1))
  expect_equal(as.numeric(table(factor(draws, U_VALUES))) / 3000, probability, tolerance = 0.05)
})

test_that("draw_cutoff returns NULL when the posterior cannot be evaluated, and stops on a missing atom", {
  tau_prior <- sampler_settings()$tau_prior
  expect_null(draw_cutoff("P", U_VALUES, FREQUENCY, modifyList(STATE, list(tau = NA_real_)),
                          tau_prior, V2021))
  broken <- modifyList(STATE, list(psi_np = c(NA, rep(0.1, 8))))
  expect_error(draw_cutoff("NP", U_VALUES, FREQUENCY, broken, tau_prior, V2021), "atom probability is missing")
})

test_that("draw_tau is a Gamma(C + kappa, zeta + 1) draw", {
  tau_prior <- list(kappa = 0.5, zeta = 1e-10)
  set.seed(6)
  draw <- draw_tau(4, tau_prior)
  set.seed(6)
  expect_equal(draw, rgamma(1, 4.5, 1 + 1e-10))
  set.seed(6)
  draws <- vapply(1:4000, function(i) draw_tau(4, tau_prior), numeric(1))
  expect_equal(mean(draws), 4.5, tolerance = 0.03)
  # without a cutoff: the prior shape
  set.seed(6)
  no_cutoff <- draw_tau(NULL, tau_prior)
  set.seed(6)
  expect_equal(no_cutoff, rgamma(1, 0.5, 1 + 1e-10))
})

test_that("draw_groups puts every count up to the cutoff in the null group", {
  set.seed(1)
  groups <- draw_groups("P", U_VALUES, X_SORTED, FREQUENCY, STATE, C = 8, V2021)
  expect_named(groups, c("K", "N0", "N1", "N2", "N3", "n00", "n0", "n1"))
  expect_equal(groups$K, 8)
  expect_equal(groups$N0, 50)
  expect_equal(groups$N1, 0)
  expect_equal(groups$N2, sum(X_SORTED))
  expect_equal(groups$N3, 0)
  expect_equal(groups$n00, 20)
  expect_equal(groups$n0, FREQUENCY)
  expect_equal(groups$n1, rep(0, 9))
})

test_that("draw_groups splits the counts above the cutoff by their local fdr", {
  # with almost no null mass the counts above the cutoff are all non-null
  state <- modifyList(STATE, list(pi0 = 1e-12))
  groups <- draw_groups("P", U_VALUES, X_SORTED, FREQUENCY, state, C = 3, V2021)
  above <- U_VALUES > 3
  expect_equal(groups$N0, sum(FREQUENCY[!above]))
  expect_equal(groups$N1, sum(FREQUENCY[above]))
  expect_equal(groups$N3, sum((U_VALUES * FREQUENCY)[above]))
  expect_equal(groups$n1, FREQUENCY * above)
  # in general the two groups add up to the data, and the draw follows the seed
  set.seed(15)
  groups <- draw_groups("NP", U_VALUES, X_SORTED, FREQUENCY, STATE, C = 3, V2021)
  expect_equal(groups$n0 + groups$n1, FREQUENCY)
  expect_equal(groups$n0[!above], FREQUENCY[!above])
  expect_equal(c(groups$N0, groups$N2), c(sum(groups$n0), sum(U_VALUES * groups$n0)))
  set.seed(15)
  expect_identical(draw_groups("NP", U_VALUES, X_SORTED, FREQUENCY, STATE, C = 3, V2021), groups)
})

test_that("draw_groups keeps a non-null count with probability 1 - fdr", {
  fdr <- mixture_fdr("NP", U_VALUES, STATE, 3, V2021)$fdr
  set.seed(23)
  null_share <- rowMeans(replicate(1500, {
    draw_groups("NP", U_VALUES, X_SORTED, FREQUENCY, STATE, C = 3, V2021)$n0 / FREQUENCY
  }))
  expect_equal(null_share, fdr, tolerance = 0.08)
})

test_that("draw_groups without a cutoff uses C = 1 or the largest count, by version", {
  set.seed(2)
  one <- draw_groups("P", U_VALUES, X_SORTED, FREQUENCY, STATE, NULL, V2021)
  set.seed(2)
  expect_equal(one, draw_groups("P", U_VALUES, X_SORTED, FREQUENCY, STATE, 1, V2021))
  # 2018: the largest count, so every count is null
  largest <- draw_groups("P", U_VALUES, X_SORTED, FREQUENCY, STATE, NULL, V2018)
  expect_equal(c(largest$N0, largest$N1), c(50, 0))
})

test_that("draw_xi is a Gamma draw with the shape and rate of the shifted Poisson", {
  groups <- list(K = 8, N1 = 10, N3 = 85)
  xi_prior <- list(kappa = 0.5, zeta = 1e-10)
  # shape N3 + kappa - (C + 1) N1 + 1 = 85 + 0.5 - 40 + 1 = 46.5, rate N1 = 10
  set.seed(3)
  draw <- draw_xi(3, groups, xi_prior)
  set.seed(3)
  expect_equal(draw, rgamma(1, shape = 46.5, rate = 10 + 1e-10))
  set.seed(3)
  draws <- vapply(1:4000, function(i) draw_xi(3, groups, xi_prior), numeric(1))
  expect_equal(mean(draws), 4.65, tolerance = 0.02)
  # without a cutoff the largest count stands in: shape 85 + 0.5 - 90 + 1 is negative
  expect_true(is.na(suppressWarnings(draw_xi(NULL, groups, xi_prior))))
})

test_that("concentration_above_cutoff is small up to the cutoff", {
  expect_equal(concentration_above_cutoff(2, 6, 0.5), c(1e-6, 1e-6, 1e-6, 0.5, 0.5, 0.5))
  expect_equal(concentration_above_cutoff(2, 6, 0.5, epsilon = 0.01), c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5))
  expect_equal(concentration_above_cutoff(-1, 4, 0.5), rep(0.5, 4))
  expect_equal(concentration_above_cutoff(9, 4, 0.5), rep(1e-6, 4))
  expect_equal(concentration_above_cutoff(0, 3, 2), c(1e-6, 2, 2))
})

test_that("draw_atoms keeps pi0 f0 up to the cutoff and spreads the rest above it", {
  psi <- rep(1 / 9, 9)
  concentration <- FREQUENCY + 1
  set.seed(10)
  atoms <- draw_atoms(U_VALUES, psi, concentration, STATE, gamma = 1, C = 3, atoms = 9, V2021)
  expect_length(atoms, 9)
  expect_equal(atoms[1:4], 0.8 * dzigp(0:3, 0.4, 1.5, 0.1))
  # the mass 1 - sum(psi[1:4]) = 5 / 9 is spread by a Dirichlet draw over the five counts above
  # the cutoff and 9 - 9 + 3 + 1 = 4 atoms of the prior without observations
  set.seed(10)
  expect_equal(atoms[5:9], (5 / 9) * rdirichlet(c(concentration[5:9], rep(1, 4)))[1:5])
  expect_true(all(atoms[5:9] > 0))
  expect_lt(sum(atoms[5:9]), 5 / 9)
  set.seed(10)
  expect_identical(draw_atoms(U_VALUES, psi, concentration, STATE, 1, 3, 9, V2021), atoms)
})

test_that("draw_atoms has more atoms than distinct counts, by the rule of the version", {
  psi <- rep(1 / 9, 9)
  # a prior of 5 atoms for 9 distinct counts: 2021 takes 9 + 10 atoms, so 19 - 9 + 3 + 1 = 14
  # without observations; 2018 takes 1.5 * 8 = 12, so 7
  set.seed(10)
  atoms <- draw_atoms(U_VALUES, psi, FREQUENCY + 1, STATE, gamma = 1, C = 3, atoms = 5, V2021)
  set.seed(10)
  expect_equal(atoms[5:9], (5 / 9) * rdirichlet(c(FREQUENCY[5:9] + 1, rep(1, 14)))[1:5])
  set.seed(10)
  atoms <- draw_atoms(U_VALUES, psi, FREQUENCY + 1, STATE, gamma = 1, C = 3, atoms = 5, V2018)
  set.seed(10)
  expect_equal(atoms[5:9], (5 / 9) * rdirichlet(c(FREQUENCY[5:9] + 1, rep(1, 7)))[1:5])
})

test_that("draw_atoms without a cutoff uses C = 0 or stops, by version", {
  psi <- rep(1 / 9, 9)
  set.seed(1)
  no_cutoff <- draw_atoms(U_VALUES, psi, FREQUENCY, STATE, 1, NULL, 5, V2021)
  set.seed(1)
  expect_equal(no_cutoff, draw_atoms(U_VALUES, psi, FREQUENCY, STATE, 1, 0, 5, V2021))
  expect_error(draw_atoms(U_VALUES, psi, FREQUENCY, STATE, 1, NULL, 5, V2018), "no cutoff could be drawn")
})

test_that("null_parameters and gp_parameters map the real line to the parameter ranges", {
  expect_equal(null_parameters(c(0, 0, 0), V2018), c(0.5, 1, 0.5))
  expect_equal(null_parameters(c(log(3), log(2), 1), V2018), c(0.75, 2, 1 / (1 + exp(-1))))
  # the 2021-02-19 version has a logistic scale of 0.5 for theta
  expect_equal(null_parameters(c(0, 0, 1), V2021), c(0.5, 1, 1 / (1 + exp(-2))))
  expect_equal(gp_parameters(c(0, 0), V2018), c(1, 0.5))
  expect_equal(gp_parameters(c(log(2.5), 1), V2021), c(2.5, 1 / (1 + exp(-2))))
})

test_that("loglik_null is the ZIGP log likelihood of the null counts", {
  n0 <- c(20, 10, 6, 3, 0, 1)
  groups <- made_up_groups(n0, n1 = rep(0, 6))
  expect_equal(loglik_null(c(0.3, 1.8, 0.2), groups), sum(n0 * log(dzigp(0:5, 0.3, 1.8, 0.2))))
  expect_equal(loglik_null(c(0.5, 1, 0), groups), sum(n0 * log(0.5 * (0:5 == 0) + 0.5 * dpois(0:5, 1))))
})

test_that("loglik_gp is the generalized Poisson log likelihood of the non-null counts", {
  n1 <- c(0, 4, 6, 3, 0, 2)
  groups <- made_up_groups(n0 = rep(0, 6), n1)
  expect_equal(loglik_gp(c(2, 0.2), groups), sum(n1 * log(dzigp(0:5, 0, 2, 0.2))))
  expect_equal(loglik_gp(c(1.5, 0), groups), sum(n1 * dpois(0:5, 1.5, log = TRUE)))
})

test_that("logpost_null adds the priors and the Jacobian to the log likelihood", {
  groups <- made_up_groups(n0 = c(20, 10, 6, 3, 0, 1), n1 = rep(0, 6))
  pre <- c(-0.4, 0.3, -1)
  par <- c(plogis(-0.4), exp(0.3), plogis(-1))
  # eta uniform, lambda^(-1/2), theta uniform on (0, 1); Jacobian eta (1 - eta) lambda theta (1 - theta)
  expected <- loglik_null(par, groups) - 0.5 * log(par[2]) +
    log(par[1] * (1 - par[1]) * par[2] * par[3] * (1 - par[3]))
  expect_equal(logpost_null(pre, groups, V2018), expected)
  # 2021: theta uniform on (0, 0.5), so a theta above 0.5 has no posterior mass
  expect_equal(logpost_null(c(-0.4, 0.3, 2), groups, V2021), -Inf)
  par <- null_parameters(pre, V2021)
  expect_equal(logpost_null(pre, groups, V2021),
               loglik_null(par, groups) - 0.5 * log(par[2]) + log(2) +
                 log(par[1] * (1 - par[1]) * par[2] * 2 * par[3] * (1 - par[3])))
})

test_that("logpost_gp adds the priors and the Jacobian to the log likelihood", {
  groups <- made_up_groups(n0 = rep(0, 6), n1 = c(0, 4, 6, 3, 0, 2))
  pre <- c(0.5, -0.3)
  par <- c(exp(0.5), plogis(-0.3))
  expected <- loglik_gp(par, groups) - 0.5 * log(par[1]) + log(par[1] * par[2] * (1 - par[2]))
  expect_equal(logpost_gp(pre, groups, V2018), expected)
  par <- gp_parameters(pre, V2021)
  expect_equal(logpost_gp(pre, groups, V2021),
               loglik_gp(par, groups) - 0.5 * log(par[1]) + log(par[1] * 2 * par[2] * (1 - par[2])))
})

test_that("positive_definite leaves a positive definite matrix and repairs a singular one", {
  M <- matrix(c(2, 0.5, 0.5, 1), 2)
  expect_equal(positive_definite(M), M)
  singular <- matrix(c(1, 1, 1, 1), 2)
  expect_error(chol(singular))
  repaired <- positive_definite(singular)
  expect_equal(repaired, singular)                   # the correction is of the size of the rounding
  expect_true(all(eigen(repaired)$values > 0))
  expect_equal(crossprod(chol(repaired)), singular)
  expect_equal(positive_definite(4), matrix(4))
  expect_error(positive_definite(matrix(1, 2, 3)), "not square")
})

test_that("adapted_proposal is the factor of the covariance of the last share of the chain", {
  set.seed(14)
  history <- matrix(rnorm(60), 20, 3)
  adapt <- list(prop = 0.5)
  factor <- adapted_proposal(history, t = 12, adapt, fallback = diag(0.1, 3))
  rows <- 6:12
  expect_equal(crossprod(factor), 6 * var(history[rows, ]) / 7)
  expect_equal(factor[lower.tri(factor)], rep(0, 3))
  # a chain that has not moved has no covariance: the starting proposal is used
  still <- matrix(0, 20, 3)
  expect_equal(adapted_proposal(still, 12, adapt, fallback = diag(0.1, 3)), chol(diag(0.1, 3)))
})

test_that("bayes_local_fdr returns a local fdr in [0, 1] for every model, the same under a seed", {
  x <- SEPARATED$x
  for (model in c("NP", "SP", "P", "GP")) {
    set.seed(1)
    fit <- bayes_local_fdr(x, model, gamma = 1, draws = 150, burn = 50)
    expect_named(fit, c("values", "frequency", "fdr", "eta", "lambda", "theta", "pi0", "tau", "C",
                        "xi", "delta", "nu", "acceptance"))
    expect_equal(fit$values, sort(unique(x)), info = model)
    expect_equal(fit$frequency, as.numeric(table(x)), info = model)
    expect_length(fit$fdr, length(fit$values))
    expect_true(all(fit$fdr >= 0 & fit$fdr <= 1), info = model)
    expect_length(fit$pi0, 150)
    expect_length(fit$acceptance, 200)
    expect_true(all(fit$eta > 0 & fit$eta < 1 & fit$lambda > 0 & fit$pi0 > 0 & fit$pi0 < 1), info = model)
    expect_true(all(fit$theta > 0 & fit$theta < 1), info = model)
    expect_true(all(fit$C %in% fit$values), info = model)
    # the null counts keep a local fdr near one, the non-null counts far above them near zero
    expect_true(all(fit$fdr[fit$values <= 2] > 0.9), info = model)
    expect_true(all(fit$fdr[fit$values >= 18] < 0.05), info = model)
    set.seed(1)
    expect_identical(bayes_local_fdr(x, model, gamma = 1, draws = 150, burn = 50), fit, info = model)
  }
})

test_that("bayes_local_fdr continues the random number stream it is given", {
  x <- SEPARATED$x
  set.seed(1)
  first <- bayes_local_fdr(x, "NP", draws = 30, burn = 10)
  second <- bayes_local_fdr(x, "NP", draws = 30, burn = 10)   # no seed in between: another chain
  expect_false(identical(first$pi0, second$pi0))
  set.seed(1)
  expect_identical(bayes_local_fdr(x, "NP", draws = 30, burn = 10)$pi0, first$pi0)
})

test_that("bayes_local_fdr follows the version asked for", {
  x <- SEPARATED$x
  fits <- lapply(names(SAMPLER_VERSIONS), function(version) {
    set.seed(1)
    bayes_local_fdr(x, "SP", gamma = 0.5, draws = 60, burn = 40, version = version)
  })
  for (fit in fits) expect_true(all(fit$fdr >= 0 & fit$fdr <= 1))
  expect_false(identical(fits[[1]]$fdr, fits[[2]]$fdr))
  expect_false(identical(fits[[2]]$fdr, fits[[3]]$fdr))
  # theta has a uniform prior on (0, 0.5) in the 2021 versions
  expect_true(all(fits[[2]]$theta < 0.5 & fits[[3]]$theta < 0.5))
})

test_that("bayes_local_fdr can fix theta at zero and run without a cutoff", {
  x <- SEPARATED$x
  set.seed(1)
  zip <- bayes_local_fdr(x, "GP", draws = 40, burn = 20, zigp = FALSE)
  expect_equal(zip$theta, rep(0, 40))
  set.seed(1)
  no_cutoff <- bayes_local_fdr(x, "NP", draws = 40, burn = 20, cutoff = FALSE)
  expect_equal(no_cutoff$C, rep(-1, 40))
  expect_equal(no_cutoff$tau, rep(0, 40))
  expect_true(all(no_cutoff$fdr >= 0 & no_cutoff$fdr <= 1))
  expect_error(bayes_local_fdr(x, "XX"), "should be one of")
})

test_that("bayes_local_fdr thins the chain", {
  x <- SEPARATED$x
  set.seed(1)
  fit <- bayes_local_fdr(x, "P", draws = 20, burn = 10, thin = 3)
  expect_length(fit$pi0, 20)
  expect_length(fit$acceptance, 70)   # 20 * 3 + 10 iterations
})
