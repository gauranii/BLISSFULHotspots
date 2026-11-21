# -----------------------------------------------------------------------------
# File:    rcpp/tests/testthat/test-rcpp.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Unit tests of the Rcpp version: every C++ function of rcpp/src/
#          against the R function of the same name, on small inputs under
#          fixed seeds. The two must give identical() values and leave the
#          random number generator in the same state; this includes whole
#          seeded chains of bayes_local_fdr() for the three versions of the
#          sampler and the four alternatives.
# Inputs:  R/*.R, rcpp/src/, data/domains/cd00083.txt.
# Outputs: None.
# Usage:   Rscript rcpp/tests/testthat.R
# -----------------------------------------------------------------------------

source("../../../R/setup.R")
source(file.path(REPO_ROOT, "tests", "testthat", "helper-inputs.R"))

# The R functions, kept apart: they call each other, never the C++ functions.
R_VERSION <- new.env(parent = globalenv())
for (file in c("distributions.R", "bayes_fdr.R", "simulate.R")) {
  sys.source(file.path(REPO_ROOT, "R", file), envir = R_VERSION)
}

# From here on the functions of the same names are the C++ functions.
source(file.path(REPO_ROOT, "rcpp", "R", "load.R"))

#' Call the R and the C++ function of a name with the same arguments under
#' the same seed, and expect the same value (or an error from both) and the
#' same state of the random number generator afterwards.
#'
#' @param name Name of the function.
#' @param ... Its arguments.
#' @param seed Seed set before each call.
#' @return The value of the R function, invisibly.
expect_same <- function(name, ..., seed = 1) {
  run <- function(f) {
    set.seed(seed)
    value <- tryCatch(suppressWarnings(f(...)), error = function(e) "an error")
    list(value = value, state = .Random.seed)
  }
  r <- run(get(name, envir = R_VERSION))
  cpp <- run(get(name, envir = globalenv()))
  expect_identical(cpp$value, r$value, label = paste(name, "in C++"))
  expect_identical(cpp$state, r$state, label = paste("random number state after", name))
  invisible(r$value)
}

VERSIONS <- names(SAMPLER_VERSIONS)
MODELS <- c("NP", "SP", "P", "GP")

# Distinct counts 0 to 8 with their frequencies, and a state of the chain.
U_VALUES <- 0:8
FREQUENCY <- c(20, 10, 6, 3, 1, 2, 3, 3, 2)
X_SORTED <- rep(U_VALUES, FREQUENCY)
STATE <- list(eta = 0.4, lambda = 1.5, theta = 0.1, pi0 = 0.8, xi = 2, delta = 2, nu = 0.2, tau = 3,
              psi_sp = c(0.3, 0.2, 0.1, 0.1, 0.05, 0.05, 0.1, 0.06, 0.04),
              psi_np = c(0.35, 0.2, 0.1, 0.05, 0.05, 0.05, 0.1, 0.06, 0.04))
# A state whose mixture density is zero at the large counts: the local fdr
# is missing there.
EMPTY_STATE <- modifyList(STATE, list(pi0 = 1, lambda = 1e-300, theta = 0,
                                      psi_sp = c(rep(1 / 3, 3), rep(0, 6)),
                                      psi_np = c(rep(1 / 3, 3), rep(0, 6))))
GROUPS <- made_up_groups(n0 = c(20, 9, 5, 2, 1, 0, 0, 0, 0), n1 = c(0, 1, 1, 1, 0, 2, 3, 3, 2), C = 3)

SEPARATED <- simulated_data(small_setting(N = 60), seed = 21)$x
CD00083 <- read_domain("cd00083", PATHS$domains)$counts

test_that("the C++ functions are in place of the R functions", {
  expect_setequal(RCPP_FUNCTIONS,
                  c("dzigp", "rzigp", "rdirichlet", "rmultinomial", "shifted_density", "mixture_fdr",
                    "draw_cutoff", "draw_tau", "draw_groups", "draw_xi",
                    "concentration_above_cutoff", "draw_atoms", "null_parameters", "gp_parameters",
                    "loglik_null", "loglik_gp", "logpost_null", "logpost_gp", "positive_definite",
                    "adapted_proposal", "bayes_local_fdr", "rmixture"))
  for (name in RCPP_FUNCTIONS) {
    cpp <- get(name, envir = globalenv())
    expect_identical(formals(cpp), formals(get(name, envir = R_VERSION)))
    expect_true(any(grepl(".Call", deparse(body(cpp)), fixed = TRUE)))
  }
})

# ----- rcpp/src/distributions.cpp -----

test_that("dzigp() equals R", {
  for (x in list(0:12, as.numeric(0:40), c(3, 0, 0, 7), numeric(0), c(0, 150, 400))) {
    expect_same("dzigp", x, 0.4, 2, 0.3)
    expect_same("dzigp", x, 0, 0.7, 0)
    expect_same("dzigp", x, 1, 25, 1)
  }
  expect_same("dzigp", 0:5, 1.2, 2, 0.3)                  # outside the range: an error
  expect_same("dzigp", 0:5, 0.4, -1, 0.3)
  expect_same("dzigp", 0:5, 0.4, 2, 1.5)
  expect_same("dzigp", 0:5, 1.2, 2, 1.5, check = FALSE)   # outside the range, not checked
  expect_same("dzigp", 0:5, 0.4, -1, 0.3, check = FALSE)
  expect_same("dzigp", c(0, NA, 2), 0.4, 2, 0.3)
})

test_that("rzigp() equals R, draw by draw", {
  for (seed in 1:60) {
    expect_same("rzigp", 0.4, 2, 0.3, seed = seed)
    expect_same("rzigp", 0, 30, 0.9, seed = seed)
    expect_same("rzigp", 0.9, 0.5, 0, seed = seed)
    expect_same("rzigp", 0.1, 400, 0.5, seed = seed)
  }
  expect_same("rzigp", 2, 1, 0.3)
  expect_same("rzigp", 0.2, 3, -0.4, check = FALSE)
})

test_that("rdirichlet() equals R", {
  for (seed in 1:10) {
    expect_same("rdirichlet", c(0.5, 1, 2, 30), seed = seed)
    expect_same("rdirichlet", rep(0.02, 25), seed = seed)
    expect_same("rdirichlet", c(1, 0, 3), seed = seed)
    expect_same("rdirichlet", c(1, NA, 3), seed = seed)
    expect_same("rdirichlet", c(1, NA, 3), TRUE, seed = seed)
  }
  expect_same("rdirichlet", numeric(0))
})

test_that("rmultinomial() equals R", {
  for (seed in 1:10) {
    expect_same("rmultinomial", 1, c(0.2, 0.3, 0.5), seed = seed)
    expect_same("rmultinomial", 40, c(0.2, 0, 0.3, 0.5, 0), seed = seed)
    expect_same("rmultinomial", 7, c(3, 1, 1, 5), seed = seed)
    expect_same("rmultinomial", 5, c(0, 0, 1), seed = seed)
    expect_same("rmultinomial", 5, 0.4, seed = seed)
    expect_same("rmultinomial", 3, c(1e-300, 1e-12, 1e-12), seed = seed)
  }
  expect_same("rmultinomial", 3, c(0.5, NA, 0.5))
  expect_same("rmultinomial", 0, c(0.5, 0.5))
})

test_that("shifted_density() equals R", {
  for (C in c(-1, 0, 3, 8, 20)) {
    for (model in c("P", "GP", "SP")) {
      expect_same("shifted_density", model, 0:12, 2.5, 1.7, 0.2, C)
      expect_same("shifted_density", model, as.numeric(c(9, 4, 30)), 0.3, 40, 0.9, C)
    }
  }
})

# ----- rcpp/src/bayes_fdr.cpp: the steps of the sampler -----

test_that("mixture_fdr() equals R", {
  for (version in SAMPLER_VERSIONS) {
    for (model in MODELS) {
      for (C in list(-1, 0, 3, 8, 30, NULL)) {
        expect_same("mixture_fdr", model, U_VALUES, STATE, C, version)
        expect_same("mixture_fdr", model, U_VALUES, EMPTY_STATE, C, version)
        expect_same("mixture_fdr", model, U_VALUES, STATE, C, version, f0 = FREQUENCY / 100)
      }
    }
  }
  short <- modifyList(STATE, list(psi_sp = STATE$psi_sp[1:5], psi_np = STATE$psi_np[1:5]))
  for (model in MODELS) expect_same("mixture_fdr", model, U_VALUES, short, 2, SAMPLER_VERSIONS[[3]])
})

test_that("draw_cutoff() equals R", {
  for (version in SAMPLER_VERSIONS) {
    for (model in MODELS) {
      for (seed in 1:5) {
        expect_same("draw_cutoff", model, U_VALUES, FREQUENCY, STATE, list(kappa = 0.5, zeta = 1e-10),
                    version, seed = seed)
        expect_same("draw_cutoff", model, as.numeric(U_VALUES), FREQUENCY, STATE,
                    list(kappa = 2, zeta = 0.5), version, seed = seed)
      }
      expect_same("draw_cutoff", model, U_VALUES, FREQUENCY, EMPTY_STATE,
                  list(kappa = 0.5, zeta = 1e-10), version)
      expect_same("draw_cutoff", model, U_VALUES, FREQUENCY, modifyList(STATE, list(tau = NaN)),
                  list(kappa = 0.5, zeta = 1e-10), version)
      expect_same("draw_cutoff", model, U_VALUES, FREQUENCY,
                  modifyList(STATE, list(psi_np = c(NA, STATE$psi_np[-1]))),
                  list(kappa = 0.5, zeta = 1e-10), version)
    }
  }
})

test_that("draw_tau() and draw_xi() equal R", {
  for (seed in 1:5) {
    for (C in list(NULL, 0, 4, 11L)) {
      expect_same("draw_tau", C, list(kappa = 0.5, zeta = 1e-10), seed = seed)
      expect_same("draw_xi", C, GROUPS, list(kappa = 0.5, zeta = 1e-10), seed = seed)
    }
    expect_same("draw_xi", -1, GROUPS, list(kappa = 3, zeta = 2), seed = seed)
  }
})

test_that("draw_groups() equals R", {
  for (version in SAMPLER_VERSIONS) {
    for (model in MODELS) {
      for (C in list(-1, 0, 3, 8, NULL)) {
        for (seed in 1:3) {
          expect_same("draw_groups", model, U_VALUES, X_SORTED, FREQUENCY, STATE, C, version,
                      seed = seed)
        }
        # A missing local fdr: missing indicators, or uniform draws in their place.
        expect_same("draw_groups", model, U_VALUES, X_SORTED, FREQUENCY, EMPTY_STATE, C, version)
      }
      expect_same("draw_groups", model, as.numeric(U_VALUES), as.numeric(X_SORTED), FREQUENCY,
                  STATE, 3, version)
    }
  }
  gaps <- c(0, 0, 0, 1, 1, 4, 4, 9, 12)   # counts that are not all observed
  expect_same("draw_groups", "NP", c(0, 1, 4, 9, 12), gaps, c(3, 2, 2, 1, 1),
              modifyList(STATE, list(psi_sp = rep(0.2, 5), psi_np = rep(0.2, 5))), 1,
              SAMPLER_VERSIONS[[1]])
})

test_that("concentration_above_cutoff() equals R", {
  for (C in c(-1, 0, 1, 4, 8, 30)) {
    expect_same("concentration_above_cutoff", C, 9, 1.5)
    expect_same("concentration_above_cutoff", C, 9, 0.02, epsilon = 1e-2)
    expect_same("concentration_above_cutoff", C, 1, 3)
  }
})

test_that("draw_atoms() equals R", {
  concentration <- 1 + GROUPS$n0 + GROUPS$n1
  for (version in SAMPLER_VERSIONS) {
    for (C in list(-1, 0, 3, 8, 30, NULL)) {
      for (seed in 1:3) {
        expect_same("draw_atoms", U_VALUES, STATE$psi_np, concentration, STATE, 1, C, 12, version,
                    seed = seed)
        # More distinct counts than atoms: the version adds atoms.
        expect_same("draw_atoms", U_VALUES, STATE$psi_np, concentration, STATE, 0.5, C, 4.5,
                    version, seed = seed)
      }
      expect_same("draw_atoms", U_VALUES, STATE$psi_np[1:5], c(concentration[1:7], NA, NA), STATE,
                  1, C, 12, version)
    }
  }
})

test_that("the parameter transforms equal R", {
  for (version in SAMPLER_VERSIONS) {
    for (pre in list(c(0, 0, 0), c(-2.5, 0.7, 1.3), c(40, -3, -40), c(800, 710, -800))) {
      expect_same("null_parameters", pre, version)
      expect_same("gp_parameters", pre[1:2], version)
      expect_same("null_parameters", matrix(pre, ncol = 1), version)
    }
  }
})

test_that("the log likelihoods and log posteriors equal R", {
  several <- list(GROUPS, made_up_groups(c(3, 1), c(0, 2), C = 0),
                  made_up_groups(c(120, 40, 22, 9, 4, 1, 0, 0, 0, 0, 0, 0),
                                 c(0, 0, 0, 0, 0, 2, 6, 9, 14, 8, 3, 1), C = 4),
                  made_up_groups(7, 0))
  for (groups in several) {
    for (par in list(c(0.4, 1.5, 0.1), c(0.01, 7, 0), c(0.99, 0.2, 0.5), c(0, 1, 1), c(1, 2, 0.2))) {
      expect_same("loglik_null", par, groups)
      expect_same("loglik_gp", par[2:3], groups)
    }
    for (version in SAMPLER_VERSIONS) {
      for (pre in list(c(0, 0, 0), c(-2.5, 0.7, 1.3), c(0.3, 1.9, -0.2), c(40, -3, -40),
                       c(800, 710, -800), c(NA, 0, 0))) {
        expect_same("logpost_null", pre, groups, version)
        expect_same("logpost_gp", pre[2:3], groups, version)
        expect_same("logpost_null", matrix(pre, ncol = 1), groups, version)
      }
    }
  }
})

test_that("positive_definite() equals R", {
  set.seed(5)
  A <- crossprod(matrix(rnorm(12), 4, 3))
  indefinite <- matrix(c(2, 1, 0, 1, -1, 0.5, 0, 0.5, 1e-9), 3, 3)
  named <- matrix(c(2, 0.3, 0.3, 1), 2, 2, dimnames = list(c("a", "b"), c("a", "b")))
  for (M in list(A, indefinite, diag(c(0.1, 0.1, 0.1)), diag(c(0.1, 0.1)), matrix(0, 3, 3),
                 matrix(1, 3, 3), A[1:2, 1:2] * 1e-12, named, matrix(4), 4, -3,
                 matrix(c(1, 2, 3, 4), 2, 2), matrix(1:6, 2, 3), c(1, 2),
                 matrix(c(1, NA, NA, 1), 2, 2))) {
    expect_same("positive_definite", M)
  }
})

test_that("adapted_proposal() equals R", {
  set.seed(8)
  walk <- apply(matrix(rnorm(900, sd = 0.2), 300, 3), 2, cumsum)
  history <- rbind(walk, matrix(NA, 100, 3))
  adapt <- sampler_settings()$adapt
  for (t in c(100, 200, 257, 300)) {
    expect_same("adapted_proposal", history, t, adapt, diag(c(0.1, 0.1, 0.1)))
    expect_same("adapted_proposal", history[, 1:2], t, adapt, diag(c(0.1, 0.1)))
  }
  expect_same("adapted_proposal", history, 300, list(prop = 0.9), diag(3))
  expect_same("adapted_proposal", history, 1, adapt, diag(3))        # rows 0 and 1
  stuck <- matrix(c(0.3, -1, 2), 200, 3, byrow = TRUE)               # a chain that never moved
  expect_same("adapted_proposal", stuck, 200, adapt, diag(c(0.1, 0.1, 0.1)))
  expect_same("adapted_proposal", history, 350, adapt, diag(c(0.1, 0.2, 0.3)))  # missing rows
  expect_same("adapted_proposal", history, 500, adapt, diag(3))      # beyond the chain: an error
})

# ----- rcpp/src/simulate.cpp -----

test_that("rmixture() equals R", {
  for (f1 in c("binomial", "geometric", "binomial mixture", "poisson mixture")) {
    for (seed in c(1, 889192, 1441)) {
      expect_same("rmixture", 200, 0.4, 2, 0.3, 0.8, f1, 5, seed)
      expect_same("rmixture", 50, 0.4, 1.5, 0, 0.35, f1, 10, seed, geometric_prob = 0.6)
    }
  }
  expect_same("rmixture", 30, 0.4, 2, 0.3, 0.8, C = 5, seed = 3)     # the default f1
  expect_same("rmixture", 30, 0.4, 2, 0.3, 0.8, "normal", 5, 3)      # unknown f1: an error
  expect_same("rmixture", 30, 0.4, 2, 1.3, 0.8, "binomial", 5, 3)
  expect_same("rmixture", 0, 0.4, 2, 0.3, 0.8, "binomial", 5, 3)
})

# ----- Whole chains -----

test_that("seeded chains of bayes_local_fdr() equal R: versions and alternatives", {
  for (version in VERSIONS) {
    for (model in MODELS) {
      fit <- expect_same("bayes_local_fdr", SEPARATED, model, 1, draws = 300, burn = 100,
                         version = version, seed = 11)
      expect_true(is.list(fit) && all(is.finite(fit$fdr)))
      fit <- expect_same("bayes_local_fdr", CD00083, model, PAPER$real_gamma[[model]],
                         sampler_settings(lambda_init = 1, theta_init = 0.5), draws = 300,
                         burn = 100, version = version, seed = PAPER$real_seed)
      expect_true(is.list(fit) && all(is.finite(fit$fdr)))
    }
  }
})

test_that("seeded chains of bayes_local_fdr() equal R: the switches and defaults", {
  for (version in VERSIONS) {
    for (model in MODELS) {
      for (switches in list(c(FALSE, TRUE), c(TRUE, FALSE), c(FALSE, FALSE))) {
        expect_same("bayes_local_fdr", SEPARATED, model, 0.5, draws = 150, burn = 60,
                    zigp = switches[1], cutoff = switches[2], version = version, seed = 4)
      }
      expect_same("bayes_local_fdr", CD00083, model, 0.02, draws = 90, burn = 45, thin = 3,
                  version = version, seed = 2)
    }
  }
  expect_same("bayes_local_fdr", SEPARATED, draws = 120, burn = 30)   # model, gamma, version by default
  expect_same("bayes_local_fdr", SEPARATED, "G", draws = 120, burn = 30)
  expect_same("bayes_local_fdr", SEPARATED, "Poisson", draws = 120, burn = 30)
  expect_same("bayes_local_fdr", SEPARATED, "NP", draws = 20, burn = 30)   # no draw to average
  expect_same("bayes_local_fdr", rep(0, 20), "NP", draws = 20, burn = 5)   # one distinct count
})

test_that("seeded chains of bayes_local_fdr() equal R: more distinct counts than atoms", {
  wide <- c(rep(0:4, c(60, 30, 14, 6, 3)), 12:70, 30:45)   # 64 distinct counts
  for (version in VERSIONS) {
    for (model in MODELS) {
      expect_same("bayes_local_fdr", wide, model, 1, draws = 150, burn = 60, version = version,
                  seed = 6)
    }
  }
})

test_that("the simulation code of R/ runs on the C++ functions and equals R", {
  setting <- small_setting(N = 80)
  cpp <- bayes_replication(setting, "NP", 1, 3, SMALL_SETTINGS)
  r_code <- bayes_replication
  environment(r_code) <- R_VERSION
  sys.source(file.path(REPO_ROOT, "R", "eb_fdr.R"), envir = R_VERSION)
  expect_identical(cpp, r_code(setting, "NP", 1, 3, SMALL_SETTINGS))
  expect_identical(domain_bayes_rejections(CD00083, "GP", SMALL_SETTINGS),
                   local({
                     r_code <- domain_bayes_rejections
                     environment(r_code) <- R_VERSION
                     r_code(CD00083, "GP", SMALL_SETTINGS)
                   }))
})
