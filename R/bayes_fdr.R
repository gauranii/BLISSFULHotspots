# -----------------------------------------------------------------------------
# File:    R/bayes_fdr.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Tables 1 to 3, Figures S1 and S3.
# Purpose: The fully Bayesian local false discovery rate: a Gibbs sampler for
#          the two-group model f = pi0 f0 + (1 - pi0) f1 with a ZIGP null f0, a
#          random cutoff C below which every count is null, and one of four
#          alternatives f1: nonparametric (NP), semiparametric (SP), Poisson (P)
#          or generalized Poisson (GP).
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R")
#          fit <- bayes_local_fdr(x, model = "NP", gamma = 1)
#          fit$values[fit$fdr <= 0.05]   # counts declared non-null
# Source:  simFDR_fcns.R (2021 simulation runs: fdr, draw_C, draw_tau, draw_z,
#          draw_pi0, draw_xi, spec_psi, tx_*, loglik_*, logpost_*, pd_matrix,
#          gibbs); fcns1_dist.R to fcns4_gibbs.R (2018 analysis of the protein
#          domains: the same functions, and spec_beta).
# -----------------------------------------------------------------------------

# Alternatives f1, in the order of the model numbers 1 to 4 of the original code.
BAYES_MODELS <- c("P", "GP", "SP", "NP")

# The sampler was revised between the analysis of the protein domains (2018)
# and the simulation runs (February 2021). Each entry lists what one version
# does where the versions differ, so that every run can be repeated as it was.
#   check_range      stop when a ZIGP parameter leaves its range
#   replace_missing  replace a missing Dirichlet concentration by 1e-5 and a
#                    missing null probability by a uniform draw
#   no_cutoff        cutoff used where none could be drawn (see draw_cutoff()):
#                    "one" or "largest count"
#   extra_atoms      atoms of the Dirichlet prior when the data have more
#                    distinct counts than the prior has atoms
#   logit_scale      scale of the logistic transform of theta and nu
#   theta_max        upper end of the uniform prior of theta
#   sp_prior         concentration of the semiparametric prior: "flat" (gamma
#                    for every count) or "above cutoff" (0.01 up to the cutoff)
#   fdr_from         first stored draw that enters the posterior mean of the fdr
SAMPLER_VERSIONS <- list(
  "2018-03-25" = list(check_range = FALSE, replace_missing = FALSE, no_cutoff = "largest count",
                      extra_atoms = function(U, largest) 1.5 * largest, logit_scale = 1,
                      theta_max = 1, sp_prior = "above cutoff", fdr_from = function(burn) 1),
  "2021-02-07" = list(check_range = TRUE, replace_missing = FALSE, no_cutoff = "one",
                      extra_atoms = function(U, largest) U + 10, logit_scale = 1,
                      theta_max = 0.5, sp_prior = "flat", fdr_from = function(burn) burn + 1),
  "2021-02-19" = list(check_range = TRUE, replace_missing = TRUE, no_cutoff = "one",
                      extra_atoms = function(U, largest) U + 10, logit_scale = 0.5,
                      theta_max = 0.5, sp_prior = "flat", fdr_from = function(burn) burn + 1)
)

#' Starting values, priors and proposal settings of the sampler.
#'
#' @param lambda_init,theta_init Starting values of the ZIGP rate and
#'   dispersion (2 and 0.25 in the simulations of Tables 1 and 2).
#' @return List with `init` (starting values; `beta` are the starting atom
#'   probabilities and `atoms` their number), `pi0_prior` (Beta), `tau_prior`
#'   and `xi_prior` (Gamma shape `kappa` and rate `zeta`), `proposal` (starting
#'   points and covariances of the two random walks) and `adapt` (the proposal
#'   covariance is re-estimated every `thin` iterations between the fractions
#'   `lower` and `upper` of the run, from the last `prop` share of the chain).
sampler_settings <- function(lambda_init = 2, theta_init = 0.25) {
  atoms <- 50
  list(init = list(eta = 0.5, lambda = lambda_init, theta = theta_init, pi0 = 0.5, xi = 1,
                   delta = 1, nu = 0.5, tau = 1, beta = rep(1 / atoms, atoms), atoms = atoms),
       pi0_prior = list(alpha = 1, beta = 1),
       tau_prior = list(kappa = 0.5, zeta = 1e-10),
       xi_prior = list(kappa = 0.5, zeta = 1e-10),
       proposal = list(null_start = c(0, 0, 0), gp_start = c(0, 0),
                       null_var = diag(c(0.1, 0.1, 0.1)), gp_var = diag(c(0.1, 0.1))),
       adapt = list(lower = 0.1, upper = 0.75, thin = 100, prop = 0.5))
}

# ----- Mixture density and local fdr -----

#' Mixture density and local false discovery rate at the distinct counts.
#'
#' f = pi0 f0 up to the cutoff; above it f = pi0 f0 + (1 - pi0) f1 (P, GP, SP)
#' or the atom probabilities psi (NP). fdr = min(1, pi0 f0 / f).
#'
#' @param model One of `BAYES_MODELS`.
#' @param u Distinct counts, increasing.
#' @param s Current state: list with eta, lambda, theta, pi0, xi, delta, nu,
#'   psi_sp and psi_np.
#' @param C Cutoff, or NULL when none could be drawn.
#' @param version Entry of `SAMPLER_VERSIONS`.
#' @param f0 pi0 f0(u), if already computed.
#' @return List with `f0` (pi0 f0), `f` and `fdr`.
mixture_fdr <- function(model, u, s, C, version,
                        f0 = s$pi0 * dzigp(u, s$eta, s$lambda, s$theta, version$check_range)) {
  if (is.null(C)) {
    if (version$no_cutoff != "one") stop("no cutoff could be drawn")
    C <- 1
  }
  f1 <- numeric(length(u))
  if (model %in% c("P", "GP")) {
    f1 <- (1 - s$pi0) * shifted_density(model, u, s$xi, s$delta, s$nu, C)
  } else if (model == "SP") {
    f1 <- (1 - s$pi0) * s$psi_sp
  }
  f <- numeric(length(u))
  below <- which(u <= C)
  above <- which(u > C)
  f[below] <- f0[below]
  if (model != "NP") {
    f[above] <- f0[above] + f1[above]
  } else {
    f[above] <- s$psi_np[above]
  }
  list(f0 = f0, f = f, fdr = pmin(1, f0 / f))
}

# ----- Full conditional draws -----

#' Draw the cutoff C from its discrete full conditional over the distinct counts.
#'
#' C has a Poisson(tau) prior and tau a Gamma prior. If the log posterior
#' cannot be evaluated no cutoff is drawn and NULL is returned.
#'
#' @param model,u,s,version As in `mixture_fdr()`.
#' @param frequency Number of observations at each distinct count.
#' @param tau_prior Gamma prior of tau (`kappa`, `zeta`).
#' @return The drawn cutoff, or NULL.
draw_cutoff <- function(model, u, frequency, s, tau_prior, version) {
  kappa <- tau_prior$kappa
  zeta <- tau_prior$zeta
  U <- length(u)
  if (anyNA(s$psi_sp) || anyNA(s$psi_np)) stop("an atom probability is missing")
  loglik <- numeric(U)
  log_poisson <- u * log(s$tau) - lgamma(u + 1) - s$tau
  log_gamma <- kappa * log(zeta) - lgamma(kappa) + (kappa - 1) * log(s$tau) - (s$tau * zeta)
  f0 <- s$pi0 * dzigp(u, s$eta, s$lambda, s$theta, version$check_range)

  for (k in 1:U) {
    log_f <- log(mixture_fdr(model, u, s, u[k], version, f0)$f)
    if (!any(is.nan(log_f)) || !any(is.na(log_f))) loglik[k] <- sum(frequency * log_f)
  }

  logpost <- loglik + log_poisson + log_gamma
  if (any(is.na(logpost))) return(NULL)
  q <- exp(logpost - max(logpost)) / sum(exp(logpost - max(logpost)))
  if (any(is.na(q))) return(NULL)
  u[which.max(rmultinomial(1, q))]
}

#' Draw tau, the mean of the Poisson prior of the cutoff.
#'
#' @param C Cutoff, or NULL.
#' @param tau_prior Gamma prior of tau.
#' @return One draw.
draw_tau <- function(C, tau_prior) {
  shape <- if (is.null(C)) tau_prior$kappa else C + tau_prior$kappa
  rgamma(1, shape, tau_prior$zeta + 1)
}

#' Draw the null indicators and count the observations in each group.
#'
#' An observation is null with probability 1 up to the cutoff and with its
#' local fdr above it.
#'
#' @param model,u,s,C,version As in `mixture_fdr()`.
#' @param x_sorted The observations, sorted.
#' @param frequency Number of observations at each distinct count.
#' @return List: `N0` and `N1` (numbers of null and of non-null observations,
#'   the latter above the cutoff), `N2` and `N3` (their sums), `n00` (null
#'   zeros), `n0` and `n1` (null and non-null observations per distinct
#'   count) and `K` (largest count).
draw_groups <- function(model, u, x_sorted, frequency, s, C, version) {
  K <- x_sorted[length(x_sorted)]
  N <- length(x_sorted)
  values <- 0:K
  if (is.null(C)) C <- if (version$no_cutoff == "one") 1 else K

  at_most_cutoff <- as.numeric(x_sorted <= C)
  fdr <- mixture_fdr(model, u, s, C, version)$fdr
  p <- pmax(at_most_cutoff, rep(fdr, times = frequency))
  if (version$replace_missing) {
    missing <- which(is.na(p))
    p[missing] <- runif(length(missing))
  }
  z <- rbinom(N, 1, p)

  if (anyNA(z)) {
    n0 <- vapply(values, function(v) sum(z * as.numeric(x_sorted == v)), numeric(1))
    n1 <- vapply(values, function(v) sum((1 - z) * as.numeric(x_sorted == v)), numeric(1))
  } else {
    n0 <- as.numeric(tabulate(x_sorted[z == 1] + 1, K + 1))
    n1 <- as.numeric(tabulate(x_sorted[z == 0] + 1, K + 1))
  }
  observed <- which(tabulate(x_sorted + 1, K + 1) > 0)
  list(K = K, N0 = sum(n0), N1 = sum(n1 * as.numeric(values > C)),
       N2 = sum(values * n0), N3 = sum(values * n1 * as.numeric(values > C)),
       n00 = n0[1], n0 = n0[observed], n1 = n1[observed])
}

#' Draw the Poisson mean xi of the shifted Poisson alternative.
#'
#' @param C Cutoff, or NULL.
#' @param groups Output of `draw_groups()`.
#' @param xi_prior Gamma prior of xi.
#' @return One draw.
draw_xi <- function(C, groups, xi_prior) {
  if (is.null(C)) C <- groups$K
  rgamma(1, shape = groups$N3 + xi_prior$kappa - (C + 1) * groups$N1 + 1,
         rate = groups$N1 + xi_prior$zeta)
}

#' Concentration of the semiparametric prior that is small up to the cutoff.
#'
#' @param C Cutoff (-1 when no cutoff is used).
#' @param P Number of atoms.
#' @param gamma Concentration above the cutoff.
#' @param epsilon Concentration up to the cutoff.
#' @return Vector of P concentrations.
concentration_above_cutoff <- function(C, P, gamma, epsilon = 1e-6) {
  if (C == -1) return(rep(gamma, P)[1:P])
  n <- length(seq(1, C + 1, 1))
  c(rep(epsilon, n), if (n > P) c() else rep(gamma, P - n))[1:P]
}

#' Draw the atom probabilities of the nonparametric part of the mixture.
#'
#' Up to the cutoff the atoms are pi0 f0; above it they are a Dirichlet draw
#' scaled to the remaining probability.
#'
#' @param u,s,C,version As in `mixture_fdr()`.
#' @param psi Current atom probabilities.
#' @param concentration Posterior Dirichlet concentration per distinct count.
#' @param gamma Prior concentration of an atom without observations.
#' @param atoms Number of atoms of the prior.
#' @return Atom probabilities at the distinct counts.
draw_atoms <- function(u, psi, concentration, s, gamma, C, atoms, version) {
  if (is.null(C)) {
    if (version$no_cutoff != "one") stop("no cutoff could be drawn")
    C <- 0
  }
  U <- length(u)
  below <- which(u <= C)
  above <- which(u > C)
  null_mass <- min(1 - 1e-10, sum(psi[below]))
  if (atoms < U) atoms <- version$extra_atoms(U, max(u))
  concentration <- c(concentration[above], rep(gamma, max(0, atoms - U + C + 1)))
  psi_null <- s$pi0 * dzigp(u[below], s$eta, s$lambda, s$theta, version$check_range)
  psi_alternative <- (1 - null_mass) * rdirichlet(concentration, version$replace_missing)
  c(psi_null, psi_alternative)[seq_len(U)]
}

# ----- Random walk Metropolis-Hastings steps -----

#' ZIGP parameters from their unconstrained scale.
#' @param pre Unconstrained (eta, lambda, theta).
#' @param version Entry of `SAMPLER_VERSIONS`.
#' @return c(eta, lambda, theta).
null_parameters <- function(pre, version) {
  c(plogis(pre[1]), exp(pre[2]), plogis(pre[3], 0, version$logit_scale))
}

#' Generalized Poisson parameters of f1 from their unconstrained scale.
#' @param pre Unconstrained (delta, nu).
#' @param version Entry of `SAMPLER_VERSIONS`.
#' @return c(delta, nu).
gp_parameters <- function(pre, version) {
  c(exp(pre[1]), plogis(pre[2], 0, version$logit_scale))
}

#' Log likelihood of the ZIGP parameters given the null observations.
#' @param par c(eta, lambda, theta).
#' @param groups Output of `draw_groups()`.
#' @return The log likelihood.
loglik_null <- function(par, groups) {
  eta <- par[1]
  lambda <- par[2]
  theta <- par[3]
  n0 <- groups$n0
  last <- length(n0)
  y <- seq(1, last - 1, 1)
  groups$n00 * log(eta + (1 - eta) * exp(-lambda)) + (groups$N0 - groups$n00) * log(1 - eta) +
    (groups$N0 - groups$n00) * (log(lambda) - lambda) - groups$N2 * theta +
    sum(n0[2:last] * (y - 1) * log(lambda + theta * y)) -
    sum(n0[2:last] * lgamma(y + 1))
}

#' Log likelihood of the generalized Poisson f1 given the non-null observations.
#' @param par c(delta, nu).
#' @param groups Output of `draw_groups()`.
#' @return The log likelihood.
loglik_gp <- function(par, groups) {
  delta <- par[1]
  nu <- par[2]
  n1 <- groups$n1
  last <- length(n1)
  y <- seq(1, last - 1, 1)
  groups$N1 * log(delta) - delta * groups$N1 - nu * groups$N3 +
    sum(n1[2:last] * (y - 1) * log(delta + nu * y)) -
    sum(n1[2:last] * lgamma(y + 1))
}

#' Log posterior of the unconstrained ZIGP parameters.
#'
#' Priors: eta uniform on (0, 1), lambda proportional to lambda^(-1/2), theta
#' uniform on (0, theta_max).
#'
#' @param pre Unconstrained (eta, lambda, theta).
#' @param groups Output of `draw_groups()`.
#' @param version Entry of `SAMPLER_VERSIONS`.
#' @return The log posterior, with the Jacobian of the transform.
logpost_null <- function(pre, groups, version) {
  par <- null_parameters(pre, version)
  prior <- dunif(par[1], log = TRUE) + -0.5 * log(par[2]) +
    dunif(par[3], 0, version$theta_max, log = TRUE)
  jacobian <- c(dlogis(pre[1]), exp(pre[2]), dlogis(pre[3], 0, version$logit_scale))
  loglik_null(par, groups) + prior + log(abs(prod(jacobian)))
}

#' Log posterior of the unconstrained generalized Poisson parameters of f1.
#'
#' Priors: delta proportional to delta^(-1/2), nu uniform on (0, 1).
#'
#' @param pre Unconstrained (delta, nu).
#' @param groups Output of `draw_groups()`.
#' @param version Entry of `SAMPLER_VERSIONS`.
#' @return The log posterior, with the Jacobian of the transform.
logpost_gp <- function(pre, groups, version) {
  par <- gp_parameters(pre, version)
  prior <- -0.5 * log(par[1]) + dunif(par[2], log = TRUE)
  jacobian <- c(exp(pre[1]), dlogis(pre[2], 0, version$logit_scale))
  loglik_gp(par, groups) + prior + log(abs(prod(jacobian)))
}

#' Nearest positive definite matrix, by lifting the small eigenvalues.
#' @param M Symmetric matrix.
#' @return `M` plus the correction.
positive_definite <- function(M) {
  if (!is.matrix(M)) M <- as.matrix(M)
  d <- dim(M)[1]
  if (dim(M)[2] != d) stop("Input matrix is not square!")
  e <- eigen(M)
  tolerance <- d * max(abs(e$values)) * .Machine$double.eps
  lift <- pmax(0, 2 * tolerance - e$values)
  M + e$vectors %*% diag(lift, d) %*% t(e$vectors)
}

#' Cholesky factor of the adapted proposal covariance.
#'
#' @param history Matrix of the chain so far, one row per iteration.
#' @param t Current iteration.
#' @param adapt Adaptation settings (`prop`).
#' @param fallback Starting proposal covariance, used if the factor fails.
#' @return Upper triangular factor.
adapted_proposal <- function(history, t, adapt, fallback) {
  rows <- floor(t * adapt$prop):t
  sigma <- (length(rows) - 1) * var(history[rows, ]) / length(rows)
  tryCatch(chol(positive_definite(sigma)),
           error = function(e) chol(positive_definite(fallback)))
}

# ----- The sampler -----

#' Bayesian local false discovery rate of count data.
#'
#' Runs the Gibbs sampler and returns the posterior mean of the local fdr at
#' every distinct count. It uses the random number stream as it finds it; set
#' the seed before the call.
#'
#' @param x Counts (non-negative integers), one per position.
#' @param model Alternative f1: "NP", "SP", "P" or "GP".
#' @param gamma Concentration of the Dirichlet prior (NP and SP).
#' @param settings Output of `sampler_settings()`.
#' @param draws Stored iterations (after burn-in, every `thin`).
#' @param burn Burn-in iterations.
#' @param thin Thinning.
#' @param zigp If FALSE theta is fixed at zero (zero-inflated Poisson null).
#' @param cutoff If FALSE no cutoff is used (C = -1).
#' @param version Name of the entry of `SAMPLER_VERSIONS` to follow.
#' @return List: `values` (distinct counts), `frequency` (observations at each),
#'   `fdr` (posterior mean local fdr at each), the stored draws `eta`,
#'   `lambda`, `theta`, `pi0`, `tau`, `C`, `xi`, `delta`, `nu`, and
#'   `acceptance` (running acceptance percentage of the ZIGP random walk).
bayes_local_fdr <- function(x, model = c("NP", "SP", "P", "GP"), gamma = 1,
                            settings = sampler_settings(), draws = 5000, burn = 1000, thin = 1,
                            zigp = TRUE, cutoff = TRUE, version = "2021-02-19") {
  model <- match.arg(model)
  version <- SAMPLER_VERSIONS[[version]]
  init <- settings$init
  adapt <- settings$adapt
  proposal <- settings$proposal

  iterations <- draws * thin + burn
  x_sorted <- sort(x)
  u <- sort(unique(x))
  U <- length(u)
  K <- max(x)
  frequency <- c(sum(x == 0), tabulate(x, K))
  frequency <- frequency[frequency > 0]

  kept <- list(eta = numeric(draws), lambda = numeric(draws), theta = numeric(draws),
               pi0 = numeric(draws), tau = numeric(draws), C = numeric(draws),
               xi = numeric(draws), delta = numeric(draws), nu = numeric(draws))
  acceptance <- numeric(iterations)
  null_history <- matrix(nrow = iterations, ncol = 3)
  gp_history <- matrix(nrow = iterations, ncol = 2)
  fdr_draws <- matrix(nrow = draws, ncol = U)

  # State of the chain.
  s <- list(eta = init$eta, lambda = init$lambda, theta = if (zigp) init$theta else 0,
            pi0 = init$pi0, xi = init$xi, delta = init$delta, nu = init$nu, tau = init$tau,
            psi_sp = init$beta[seq_len(U)], psi_np = init$beta[seq_len(U)])
  atoms <- min(1.5 * K, init$atoms)
  if (any(is.na(s$psi_sp))) s$psi_sp <- rdirichlet(rep(0.5, atoms), version$replace_missing)
  if (any(is.na(s$psi_np))) s$psi_np <- rdirichlet(rep(0.5, atoms), version$replace_missing)

  null_pre <- proposal$null_start
  gp_pre <- proposal$gp_start
  null_half <- chol(proposal$null_var)
  gp_half <- chol(proposal$gp_var)
  accepted <- 0
  h <- 1

  for (t in 1:iterations) {
    store <- (t %% thin == 0) && (t > burn)
    adapt_now <- (t > adapt$lower * iterations) && (t %% adapt$thin == 0) &&
      (t < adapt$upper * iterations)

    # Cutoff C and its prior mean tau.
    if (cutoff) {
      C <- draw_cutoff(model, u, frequency, s, settings$tau_prior, version)
      s$tau <- draw_tau(C, settings$tau_prior)
    } else {
      C <- -1
      s$tau <- 0
    }

    # Group indicators and the null proportion pi0.
    groups <- draw_groups(model, u, x_sorted, frequency, s, C, version)
    s$pi0 <- rbeta(1, groups$N0 + settings$pi0_prior$alpha, groups$N1 + settings$pi0_prior$beta)

    # Posterior concentration of the atoms.
    if (model == "SP") {
      prior <- if (version$sp_prior == "flat") rep(gamma, U)
               else concentration_above_cutoff(C, U, gamma, epsilon = 1e-2)
      concentration <- prior + groups$n1
    } else if (model == "NP") {
      concentration <- rep(gamma, U) + groups$n0 + groups$n1
    }

    # ZIGP parameters eta, lambda, theta: adaptive random walk.
    null_history[t, ] <- null_pre
    if (adapt_now) null_half <- adapted_proposal(null_history, t, adapt, proposal$null_var)
    candidate <- t(null_half) %*% rnorm(3) + null_pre
    log_ratio <- logpost_null(candidate, groups, version) - logpost_null(null_pre, groups, version)
    if (!is.na(log_ratio)) {
      if (log(runif(1)) < log_ratio) {
        null_pre <- candidate
        accepted <- accepted + 1
      }
    }
    par <- null_parameters(null_pre, version)
    s$eta <- par[1]
    s$lambda <- min(par[2], K)
    s$theta <- if (zigp) par[3] else 0
    acceptance[t] <- accepted * 100 / t

    # Parameters of a parametric f1.
    if (model == "P") {
      s$xi <- draw_xi(C, groups, settings$xi_prior)
    } else if (model == "GP") {
      gp_history[t, ] <- gp_pre
      if (adapt_now) gp_half <- adapted_proposal(gp_history, t, adapt, proposal$gp_var)
      candidate <- t(gp_half) %*% rnorm(2) + gp_pre
      log_ratio <- logpost_gp(candidate, groups, version) - logpost_gp(gp_pre, groups, version)
      if (!is.na(log_ratio)) {
        if (log(runif(1)) < log_ratio) gp_pre <- candidate
      }
      par <- gp_parameters(gp_pre, version)
      s$delta <- min(par[1], K)
      s$nu <- par[2]
    }

    # Atom probabilities of a semiparametric or nonparametric f1.
    if (model == "SP") {
      s$psi_sp <- draw_atoms(u, s$psi_sp, concentration, s, gamma, C, atoms, version)
    } else if (model == "NP") {
      s$psi_np <- draw_atoms(u, s$psi_np, concentration, s, gamma, C, atoms, version)
    }

    fdr <- mixture_fdr(model, u, s, C, version)$fdr
    if (store) {
      if (is.null(C)) stop("no cutoff could be drawn")
      kept$C[h] <- C
      for (name in c("tau", "pi0", "lambda", "eta", "theta", "xi", "delta", "nu")) {
        kept[[name]][h] <- s[[name]]
      }
      fdr_draws[h, ] <- fdr
      h <- h + 1
    }
  }

  c(list(values = u, frequency = frequency,
         fdr = colMeans(fdr_draws[version$fdr_from(burn):draws, ])),
    kept, list(acceptance = acceptance))
}
