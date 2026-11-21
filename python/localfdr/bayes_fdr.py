# -----------------------------------------------------------------------------
# File:    python/localfdr/bayes_fdr.py
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
# Usage:   from localfdr import rng
#          from localfdr.bayes_fdr import bayes_local_fdr
#          rng.set_seed(1); fit = bayes_local_fdr(x, model="NP", gamma=1)
#          fit["values"][fit["fdr"] <= 0.05]   # counts declared non-null
# Source:  R/bayes_fdr.R. A seeded chain gives the same draws as in R: the
#          random numbers come from rng.py, the Cholesky factor and the
#          eigenvalues of the adapted proposal from the LAPACK library R uses
#          (rbase.py), and every sum and product is rounded in R's order.
#          A missing value (R's NA or NaN) is NaN here, and a cutoff that
#          could not be drawn (R's NULL) is None.
# -----------------------------------------------------------------------------
import math

import numpy as np
from numba import njit

from . import rbase, rmath, rng
from .distributions import (_dzigp, _log, _rdirichlet, _rmultinomial, _shifted_density,
                            check_zigp, dzigp)
from .rmath import _dlogis, _dunif_log, _lgammafn, _plogis
from .rng import _rbinom, _unif

# Alternatives f1, in the order of the model numbers 1 to 4 of the original code.
BAYES_MODELS = ("P", "GP", "SP", "NP")

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
SAMPLER_VERSIONS = {
    "2018-03-25": {"check_range": False, "replace_missing": False, "no_cutoff": "largest count",
                   "extra_atoms": lambda U, largest: 1.5 * largest, "logit_scale": 1,
                   "theta_max": 1, "sp_prior": "above cutoff", "fdr_from": lambda burn: 1},
    "2021-02-07": {"check_range": True, "replace_missing": False, "no_cutoff": "one",
                   "extra_atoms": lambda U, largest: U + 10, "logit_scale": 1,
                   "theta_max": 0.5, "sp_prior": "flat", "fdr_from": lambda burn: burn + 1},
    "2021-02-19": {"check_range": True, "replace_missing": True, "no_cutoff": "one",
                   "extra_atoms": lambda U, largest: U + 10, "logit_scale": 0.5,
                   "theta_max": 0.5, "sp_prior": "flat", "fdr_from": lambda burn: burn + 1},
}


def sampler_settings(lambda_init=2, theta_init=0.25):
    """Starting values, priors and proposal settings of the sampler.

    Parameters: lambda_init, theta_init (starting values of the ZIGP rate
    and dispersion; 2 and 0.25 in the simulations of Tables 1 and 2).
    Returns: dict with `init` (starting values; `beta` are the starting atom
    probabilities and `atoms` their number), `pi0_prior` (Beta), `tau_prior`
    and `xi_prior` (Gamma shape `kappa` and rate `zeta`), `proposal`
    (starting points and covariances of the two random walks) and `adapt`
    (the proposal covariance is re-estimated every `thin` iterations between
    the fractions `lower` and `upper` of the run, from the last `prop` share
    of the chain).
    """
    atoms = 50
    return {"init": {"eta": 0.5, "lambda": lambda_init, "theta": theta_init, "pi0": 0.5, "xi": 1,
                     "delta": 1, "nu": 0.5, "tau": 1, "beta": np.full(atoms, 1 / atoms),
                     "atoms": atoms},
            "pi0_prior": {"alpha": 1, "beta": 1},
            "tau_prior": {"kappa": 0.5, "zeta": 1e-10},
            "xi_prior": {"kappa": 0.5, "zeta": 1e-10},
            "proposal": {"null_start": np.zeros(3), "gp_start": np.zeros(2),
                         "null_var": np.diag([0.1, 0.1, 0.1]), "gp_var": np.diag([0.1, 0.1])},
            "adapt": {"lower": 0.1, "upper": 0.75, "thin": 100, "prop": 0.5}}


def _model_code(model):
    if model not in BAYES_MODELS:
        raise ValueError(f"'model' should be one of {', '.join(BAYES_MODELS)}")
    return BAYES_MODELS.index(model)


def _vector(x):
    return np.atleast_1d(np.asarray(x, dtype=np.float64))


def _null_density(u, s, version):
    """pi0 f0(u)."""
    return s["pi0"] * dzigp(u, s["eta"], s["lambda"], s["theta"], version["check_range"])


# ----- Mixture density and local fdr -----

@njit(cache=True, error_model="numpy")
def _mixture(code, u, f0, pi0, xi, delta, nu, psi, C, f, f1):
    """Fill f with the mixture density at the distinct counts u for cutoff C
    (f1 is work space). An atom beyond the end of psi is missing; returns
    whether a missing atom was used where the rest of the sum is a number."""
    used_missing_atom = False
    if code <= 1:
        _shifted_density(code == 0, u, xi, delta, nu, C, f1)
    for i in range(len(u)):
        if u[i] <= C:
            f[i] = f0[i]
        elif u[i] > C:
            if code <= 1:
                f[i] = f0[i] + (1 - pi0) * f1[i]
            elif code == 2:
                if i < len(psi):
                    f[i] = f0[i] + (1 - pi0) * psi[i]
                else:
                    used_missing_atom = used_missing_atom or not math.isnan(f0[i])
                    f[i] = np.nan
            else:
                if i < len(psi):
                    f[i] = psi[i]
                else:
                    used_missing_atom = True
                    f[i] = np.nan
        else:
            f[i] = 0.0
    return used_missing_atom


def mixture_fdr(model, u, s, C, version, f0=None):
    """Mixture density and local false discovery rate at the distinct counts.

    f = pi0 f0 up to the cutoff; above it f = pi0 f0 + (1 - pi0) f1 (P, GP,
    SP) or the atom probabilities psi (NP). fdr = min(1, pi0 f0 / f).

    Parameters: model (one of `BAYES_MODELS`), u (distinct counts,
    increasing), s (current state: dict with eta, lambda, theta, pi0, xi,
    delta, nu, psi_sp and psi_np), C (cutoff, or None when none could be
    drawn), version (entry of `SAMPLER_VERSIONS`), f0 (pi0 f0(u), if already
    computed).
    Returns: dict with `f0` (pi0 f0), `f` and `fdr`.
    """
    code = _model_code(model)
    u = _vector(u)
    if f0 is None:
        f0 = _null_density(u, s, version)
    if C is None:
        if version["no_cutoff"] != "one":
            raise ValueError("no cutoff could be drawn")
        C = 1
    psi = _vector(s["psi_sp"] if code == 2 else s["psi_np"])
    f = np.empty(len(u))
    _mixture(code, u, f0, float(s["pi0"]), float(s["xi"]), float(s["delta"]), float(s["nu"]), psi,
             float(C), f, np.empty(len(u)))
    with np.errstate(divide="ignore", invalid="ignore"):
        fdr = rbase.pmin(1.0, f0 / f)
    return {"f0": f0, "f": f, "fdr": fdr}


# ----- Full conditional draws -----

@njit(cache=True, error_model="numpy")
def _cutoff_probabilities(code, u, frequency, f0, pi0, xi, delta, nu, psi, tau, kappa, zeta, q):
    """Fill q with the full conditional probabilities of the cutoff over the
    distinct counts; returns False when they cannot be evaluated."""
    U = len(u)
    f = np.empty(U)
    f1 = np.empty(U)
    log_gamma = kappa * _log(zeta) - _lgammafn(kappa) + (kappa - 1) * _log(tau) - (tau * zeta)
    largest = -np.inf
    for k in range(U):
        used_missing_atom = _mixture(code, u, f0, pi0, xi, delta, nu, psi, u[k], f, f1)
        not_a_number = False
        loglik = 0.0
        for i in range(U):
            log_f = _log(f[i])
            if math.isnan(log_f) and not (i >= len(psi) and code >= 2 and u[i] > u[k]
                                          and (code == 3 or not math.isnan(f0[i]))):
                not_a_number = True
            loglik += frequency[i] * log_f
        if not_a_number:
            loglik = 0.0
        elif used_missing_atom:
            loglik = np.nan
        log_poisson = u[k] * _log(tau) - _lgammafn(u[k] + 1) - tau
        q[k] = loglik + log_poisson + log_gamma
        if math.isnan(q[k]):
            return False
        if q[k] > largest:
            largest = q[k]
    total = 0.0
    for k in range(U):
        q[k] = math.exp(q[k] - largest)
        total += q[k]
    for k in range(U):
        q[k] = q[k] / total
        if math.isnan(q[k]):
            return False
    return True


def draw_cutoff(model, u, frequency, s, tau_prior, version):
    """Draw the cutoff C from its discrete full conditional over the distinct
    counts.

    C has a Poisson(tau) prior and tau a Gamma prior. If the log posterior
    cannot be evaluated no cutoff is drawn and None is returned.

    Parameters: model, u, s, version (as in `mixture_fdr()`), frequency
    (number of observations at each distinct count), tau_prior (Gamma prior
    of tau: `kappa`, `zeta`).
    Returns: the drawn cutoff, or None.
    """
    code = _model_code(model)
    u = _vector(u)
    if np.isnan(s["psi_sp"]).any() or np.isnan(s["psi_np"]).any():
        raise ValueError("an atom probability is missing")
    f0 = _null_density(u, s, version)
    psi = _vector(s["psi_sp"] if code == 2 else s["psi_np"])
    q = np.empty(len(u))
    if not _cutoff_probabilities(code, u, _vector(frequency), f0, float(s["pi0"]), float(s["xi"]),
                                 float(s["delta"]), float(s["nu"]), psi, float(s["tau"]),
                                 float(tau_prior["kappa"]), float(tau_prior["zeta"]), q):
        return None
    draw = np.empty(len(u))
    if not _rmultinomial(rng._state, 1.0, q, draw):
        raise ValueError("missing value where TRUE/FALSE needed")
    return float(u[np.argmax(draw)])


def draw_tau(C, tau_prior):
    """Draw tau, the mean of the Poisson prior of the cutoff.

    Parameters: C (cutoff, or None), tau_prior (Gamma prior of tau).
    Returns: one draw (float).
    """
    shape = tau_prior["kappa"] if C is None else C + tau_prior["kappa"]
    return float(rng.rgamma(1, shape, tau_prior["zeta"] + 1)[0])


@njit(cache=True, error_model="numpy")
def _groups(state, x_sorted, fdr, frequency, C, replace_missing, n0, n1):
    """Draw the null indicators; fill n0 and n1 (observations per count 0,
    ..., K in each group). Returns False when an indicator is missing."""
    N = len(x_sorted)
    p = np.empty(N)
    position = 0
    for k in range(len(fdr)):
        for _ in range(int(frequency[k])):
            at_most_cutoff = 1.0 if x_sorted[position] <= C else 0.0
            if math.isnan(fdr[k]):
                p[position] = np.nan
            else:
                p[position] = fdr[k] if fdr[k] > at_most_cutoff else at_most_cutoff
            position += 1
    if replace_missing:
        for i in range(N):
            if math.isnan(p[i]):
                p[i] = _unif(state)
    n0[:] = 0.0
    n1[:] = 0.0
    complete = True
    for i in range(N):
        z = _rbinom(state, 1.0, p[i])
        if math.isnan(z):
            complete = False
        elif z == 1:
            n0[int(x_sorted[i])] += 1
        else:
            n1[int(x_sorted[i])] += 1
    return complete


def draw_groups(model, u, x_sorted, frequency, s, C, version):
    """Draw the null indicators and count the observations in each group.

    An observation is null with probability 1 up to the cutoff and with its
    local fdr above it.

    Parameters: model, u, s, C, version (as in `mixture_fdr()`), x_sorted
    (the observations, sorted), frequency (number of observations at each
    distinct count).
    Returns: dict: `N0` and `N1` (numbers of null and of non-null
    observations, the latter above the cutoff), `N2` and `N3` (their sums),
    `n00` (null zeros), `n0` and `n1` (null and non-null observations per
    distinct count) and `K` (largest count).
    """
    u = _vector(u)
    x_sorted = _vector(x_sorted)
    frequency = _vector(frequency)
    K = float(x_sorted[-1])
    if C is None:
        C = 1 if version["no_cutoff"] == "one" else K
    fdr = mixture_fdr(model, u, s, C, version)["fdr"]
    if len(fdr) != len(frequency) or frequency.sum() != len(x_sorted):
        raise ValueError("invalid 'times' argument")
    n0 = np.empty(int(K) + 1)
    n1 = np.empty(int(K) + 1)
    complete = _groups(rng._state, x_sorted, fdr, frequency, float(C),
                       bool(version["replace_missing"]), n0, n1)
    if not complete:
        n0[:] = np.nan
        n1[:] = np.nan
    values = np.arange(int(K) + 1, dtype=np.float64)
    above = (values > C).astype(np.float64)
    observed = np.unique(x_sorted).astype(np.int64)
    return {"K": K, "N0": rbase.r_sum(n0), "N1": rbase.r_sum(n1 * above),
            "N2": rbase.r_sum(values * n0), "N3": rbase.r_sum(values * n1 * above),
            "n00": float(n0[0]), "n0": n0[observed], "n1": n1[observed]}


def draw_xi(C, groups, xi_prior):
    """Draw the Poisson mean xi of the shifted Poisson alternative.

    Parameters: C (cutoff, or None), groups (output of `draw_groups()`),
    xi_prior (Gamma prior of xi).
    Returns: one draw (float).
    """
    if C is None:
        C = groups["K"]
    shape = groups["N3"] + xi_prior["kappa"] - (C + 1) * groups["N1"] + 1
    return float(rng.rgamma(1, shape, groups["N1"] + xi_prior["zeta"])[0])


def concentration_above_cutoff(C, P, gamma, epsilon=1e-6):
    """Concentration of the semiparametric prior that is small up to the
    cutoff.

    Parameters: C (cutoff; -1 when no cutoff is used), P (number of atoms),
    gamma (concentration above the cutoff), epsilon (concentration up to
    the cutoff).
    Returns: array of P concentrations.
    """
    if C is None:
        raise ValueError("argument is of length zero")
    P = int(P)
    if C == -1:
        return np.full(P, float(gamma))
    n = int(C) + 1
    return np.concatenate((np.full(n, float(epsilon)), np.full(max(P - n, 0), float(gamma))))[:P]


def draw_atoms(u, psi, concentration, s, gamma, C, atoms, version):
    """Draw the atom probabilities of the nonparametric part of the mixture.

    Up to the cutoff the atoms are pi0 f0; above it they are a Dirichlet draw
    scaled to the remaining probability.

    Parameters: u, s, C, version (as in `mixture_fdr()`), psi (current atom
    probabilities), concentration (posterior Dirichlet concentration per
    distinct count), gamma (prior concentration of an atom without
    observations), atoms (number of atoms of the prior).
    Returns: atom probabilities at the distinct counts.
    """
    if C is None:
        if version["no_cutoff"] != "one":
            raise ValueError("no cutoff could be drawn")
        C = 0
    u = _vector(u)
    psi = _vector(psi)
    concentration = _vector(concentration)
    U = len(u)
    below = u <= C
    n_below = int(below.sum())  # the distinct counts are increasing
    null_mass = rbase.r_min(1 - 1e-10, rbase.r_sum(psi[:n_below]) if n_below <= len(psi) else np.nan)
    if atoms < U:
        atoms = version["extra_atoms"](U, u.max())
    extra = int(max(0, atoms - U + C + 1))
    concentration = np.concatenate((concentration[~below], np.full(extra, float(gamma))))
    psi_null = s["pi0"] * dzigp(u[below], s["eta"], s["lambda"], s["theta"], version["check_range"])
    psi_alternative = (1 - null_mass) * _rdirichlet(rng._state, concentration,
                                                    bool(version["replace_missing"]))
    drawn = np.concatenate((psi_null, psi_alternative))
    if len(drawn) < U:
        drawn = np.concatenate((drawn, np.full(U - len(drawn), np.nan)))
    return drawn[:U]


# ----- Random walk Metropolis-Hastings steps -----

def null_parameters(pre, version):
    """ZIGP parameters from their unconstrained scale.

    Parameters: pre (unconstrained eta, lambda, theta), version (entry of
    `SAMPLER_VERSIONS`).
    Returns: array (eta, lambda, theta).
    """
    pre = np.asarray(pre, dtype=np.float64).ravel()
    return np.array([rmath.plogis(pre[0]), rmath.exp(pre[1]),
                     rmath.plogis(pre[2], 0, version["logit_scale"])])


def gp_parameters(pre, version):
    """Generalized Poisson parameters of f1 from their unconstrained scale.

    Parameters: pre (unconstrained delta, nu), version.
    Returns: array (delta, nu).
    """
    pre = np.asarray(pre, dtype=np.float64).ravel()
    return np.array([rmath.exp(pre[0]), rmath.plogis(pre[1], 0, version["logit_scale"])])


@njit(cache=True, error_model="numpy")
def _loglik_null(eta, lam, theta, n00, N0, N2, n0):
    with_rate = 0.0
    with_factorial = 0.0
    for y in range(1, len(n0)):
        with_rate += n0[y] * (y - 1) * _log(lam + theta * y)
        with_factorial += n0[y] * _lgammafn(y + 1.0)
    return n00 * _log(eta + (1 - eta) * math.exp(-lam)) + (N0 - n00) * _log(1 - eta) + \
        (N0 - n00) * (_log(lam) - lam) - N2 * theta + with_rate - with_factorial


@njit(cache=True, error_model="numpy")
def _loglik_gp(delta, nu, N1, N3, n1):
    with_rate = 0.0
    with_factorial = 0.0
    for y in range(1, len(n1)):
        with_rate += n1[y] * (y - 1) * _log(delta + nu * y)
        with_factorial += n1[y] * _lgammafn(y + 1.0)
    return N1 * _log(delta) - delta * N1 - nu * N3 + with_rate - with_factorial


@njit(cache=True, error_model="numpy")
def _logpost_null(pre0, pre1, pre2, n00, N0, N2, n0, logit_scale, theta_max):
    eta = _plogis(pre0, 0.0, 1.0)
    lam = math.exp(pre1)
    theta = _plogis(pre2, 0.0, logit_scale)
    prior = _dunif_log(eta, 0.0, 1.0) + -0.5 * _log(lam) + _dunif_log(theta, 0.0, theta_max)
    jacobian = _dlogis(pre0, 0.0, 1.0) * math.exp(pre1) * _dlogis(pre2, 0.0, logit_scale)
    return _loglik_null(eta, lam, theta, n00, N0, N2, n0) + prior + _log(abs(jacobian))


@njit(cache=True, error_model="numpy")
def _logpost_gp(pre0, pre1, N1, N3, n1, logit_scale):
    delta = math.exp(pre0)
    nu = _plogis(pre1, 0.0, logit_scale)
    prior = -0.5 * _log(delta) + _dunif_log(nu, 0.0, 1.0)
    jacobian = math.exp(pre0) * _dlogis(pre1, 0.0, logit_scale)
    return _loglik_gp(delta, nu, N1, N3, n1) + prior + _log(abs(jacobian))


def _counts(groups, name):
    counts = _vector(groups[name])
    if len(counts) < 2:
        raise ValueError("wrong sign in 'by' argument")  # seq(1, 0, 1) in R
    return counts


def loglik_null(par, groups):
    """Log likelihood of the ZIGP parameters given the null observations.

    Parameters: par (eta, lambda, theta), groups (output of `draw_groups()`).
    Returns: the log likelihood (float).
    """
    return float(_loglik_null(float(par[0]), float(par[1]), float(par[2]), float(groups["n00"]),
                              float(groups["N0"]), float(groups["N2"]), _counts(groups, "n0")))


def loglik_gp(par, groups):
    """Log likelihood of the generalized Poisson f1 given the non-null
    observations.

    Parameters: par (delta, nu), groups (output of `draw_groups()`).
    Returns: the log likelihood (float).
    """
    return float(_loglik_gp(float(par[0]), float(par[1]), float(groups["N1"]),
                            float(groups["N3"]), _counts(groups, "n1")))


def logpost_null(pre, groups, version):
    """Log posterior of the unconstrained ZIGP parameters.

    Priors: eta uniform on (0, 1), lambda proportional to lambda^(-1/2),
    theta uniform on (0, theta_max).

    Parameters: pre (unconstrained eta, lambda, theta), groups (output of
    `draw_groups()`), version (entry of `SAMPLER_VERSIONS`).
    Returns: the log posterior, with the Jacobian of the transform.
    """
    pre = np.asarray(pre, dtype=np.float64).ravel()
    return float(_logpost_null(pre[0], pre[1], pre[2], float(groups["n00"]), float(groups["N0"]),
                               float(groups["N2"]), _counts(groups, "n0"),
                               float(version["logit_scale"]), float(version["theta_max"])))


def logpost_gp(pre, groups, version):
    """Log posterior of the unconstrained generalized Poisson parameters of f1.

    Priors: delta proportional to delta^(-1/2), nu uniform on (0, 1).

    Parameters: pre (unconstrained delta, nu), groups, version.
    Returns: the log posterior, with the Jacobian of the transform.
    """
    pre = np.asarray(pre, dtype=np.float64).ravel()
    return float(_logpost_gp(pre[0], pre[1], float(groups["N1"]), float(groups["N3"]),
                             _counts(groups, "n1"), float(version["logit_scale"])))


def positive_definite(M):
    """Nearest positive definite matrix, by lifting the small eigenvalues.

    Parameters: M (symmetric matrix). Returns: `M` plus the correction.
    """
    M = rbase.as_matrix(M)
    d = M.shape[0]
    if M.shape[1] != d:
        raise ValueError("Input matrix is not square!")
    e = rbase.eigen_symmetric(M)
    tolerance = d * np.max(np.abs(e["values"])) * np.finfo(float).eps
    lift = np.maximum(0.0, 2 * tolerance - e["values"])
    return M + rbase.matprod(rbase.matprod(e["vectors"], np.diag(lift)), e["vectors"].T)


def adapted_proposal(history, t, adapt, fallback):
    """Cholesky factor of the adapted proposal covariance.

    Parameters: history (matrix of the chain so far, one row per iteration),
    t (current iteration, counted from 1), adapt (adaptation settings,
    `prop`), fallback (starting proposal covariance, used if the factor
    fails).
    Returns: upper triangular factor.
    """
    first = math.floor(t * adapt["prop"])
    n_rows = t - first + 1
    sigma = (n_rows - 1) * rbase.r_var(np.asarray(history)[max(first, 1) - 1:t]) / n_rows
    try:
        return rbase.chol(positive_definite(sigma))
    except ValueError:
        return rbase.chol(positive_definite(fallback))


def _proposal(half, start):
    """t(half) %*% rnorm(length(start)) + start."""
    step = rbase.matprod(np.asarray(half).T, rng.rnorm(len(start)))[:, 0]
    return step + start


# ----- The sampler -----

def bayes_local_fdr(x, model="NP", gamma=1, settings=None, draws=5000, burn=1000, thin=1,
                    zigp=True, cutoff=True, version="2021-02-19"):
    """Bayesian local false discovery rate of count data.

    Runs the Gibbs sampler and returns the posterior mean of the local fdr at
    every distinct count. It uses the random number stream as it finds it;
    set the seed (`rng.set_seed()`) before the call.

    Parameters: x (counts, non-negative integers, one per position), model
    (alternative f1: "NP", "SP", "P" or "GP"), gamma (concentration of the
    Dirichlet prior, NP and SP), settings (output of `sampler_settings()`),
    draws (stored iterations, after burn-in, every `thin`), burn (burn-in
    iterations), thin (thinning), zigp (if False theta is fixed at zero:
    zero-inflated Poisson null), cutoff (if False no cutoff is used, C = -1),
    version (name of the entry of `SAMPLER_VERSIONS` to follow).
    Returns: dict: `values` (distinct counts), `frequency` (observations at
    each), `fdr` (posterior mean local fdr at each), the stored draws `eta`,
    `lambda`, `theta`, `pi0`, `tau`, `C`, `xi`, `delta`, `nu`, and
    `acceptance` (running acceptance percentage of the ZIGP random walk).
    """
    code = _model_code(model)
    version = SAMPLER_VERSIONS[version]
    settings = sampler_settings() if settings is None else settings
    init = settings["init"]
    adapt = settings["adapt"]
    proposal = settings["proposal"]
    logit_scale, theta_max = float(version["logit_scale"]), float(version["theta_max"])

    iterations = draws * thin + burn
    x = _vector(x)
    x_sorted = np.sort(x)
    u = np.unique(x)
    U = len(u)
    K = float(x.max())
    frequency = rbase.count_frequencies_from_zero(x, K)
    frequency = frequency[frequency > 0]

    names = ("eta", "lambda", "theta", "pi0", "tau", "C", "xi", "delta", "nu")
    kept = {name: np.zeros(draws) for name in names}
    acceptance = np.zeros(iterations)
    null_history = np.full((iterations, 3), np.nan)
    gp_history = np.full((iterations, 2), np.nan)
    fdr_draws = np.full((draws, U), np.nan)

    # State of the chain.
    beta = _vector(init["beta"])
    start = np.concatenate((beta[:U], np.full(max(U - len(beta), 0), np.nan)))
    s = {"eta": init["eta"], "lambda": init["lambda"], "theta": init["theta"] if zigp else 0,
         "pi0": init["pi0"], "xi": init["xi"], "delta": init["delta"], "nu": init["nu"],
         "tau": init["tau"], "psi_sp": start, "psi_np": start.copy()}
    atoms = min(1.5 * K, init["atoms"])
    for name in ("psi_sp", "psi_np"):
        if np.isnan(s[name]).any():
            s[name] = _rdirichlet(rng._state, np.full(int(atoms), 0.5),
                                  bool(version["replace_missing"]))

    null_pre = _vector(proposal["null_start"])
    gp_pre = _vector(proposal["gp_start"])
    null_half = rbase.chol(proposal["null_var"])
    gp_half = rbase.chol(proposal["gp_var"])
    accepted = 0
    h = 0
    concentration = None

    for t in range(1, iterations + 1):
        store = (t % thin == 0) and (t > burn)
        adapt_now = (t > adapt["lower"] * iterations) and (t % adapt["thin"] == 0) and \
            (t < adapt["upper"] * iterations)

        # Cutoff C and its prior mean tau.
        if cutoff:
            C = draw_cutoff(model, u, frequency, s, settings["tau_prior"], version)
            s["tau"] = draw_tau(C, settings["tau_prior"])
        else:
            C = -1
            s["tau"] = 0

        # Group indicators and the null proportion pi0.
        groups = draw_groups(model, u, x_sorted, frequency, s, C, version)
        s["pi0"] = float(rng.rbeta(1, groups["N0"] + settings["pi0_prior"]["alpha"],
                                   groups["N1"] + settings["pi0_prior"]["beta"])[0])
        n0, n1 = _counts(groups, "n0"), _counts(groups, "n1")

        # Posterior concentration of the atoms.
        if code == 2:
            prior = np.full(U, float(gamma)) if version["sp_prior"] == "flat" else \
                concentration_above_cutoff(C, U, gamma, epsilon=1e-2)
            concentration = prior + n1
        elif code == 3:
            concentration = np.full(U, float(gamma)) + n0 + n1

        # ZIGP parameters eta, lambda, theta: adaptive random walk.
        null_history[t - 1] = null_pre
        if adapt_now:
            null_half = adapted_proposal(null_history, t, adapt, proposal["null_var"])
        candidate = _proposal(null_half, null_pre)
        log_ratio = _logpost_null(candidate[0], candidate[1], candidate[2], groups["n00"],
                                  groups["N0"], groups["N2"], n0, logit_scale, theta_max) - \
            _logpost_null(null_pre[0], null_pre[1], null_pre[2], groups["n00"], groups["N0"],
                          groups["N2"], n0, logit_scale, theta_max)
        if not math.isnan(log_ratio):
            if rmath.log(rng.runif(1)[0]) < log_ratio:
                null_pre = candidate
                accepted += 1
        par = null_parameters(null_pre, version)
        s["eta"] = float(par[0])
        s["lambda"] = rbase.r_min(par[1], K)
        s["theta"] = float(par[2]) if zigp else 0
        acceptance[t - 1] = accepted * 100 / t

        # Parameters of a parametric f1.
        if code == 0:
            s["xi"] = draw_xi(C, groups, settings["xi_prior"])
        elif code == 1:
            gp_history[t - 1] = gp_pre
            if adapt_now:
                gp_half = adapted_proposal(gp_history, t, adapt, proposal["gp_var"])
            candidate = _proposal(gp_half, gp_pre)
            log_ratio = _logpost_gp(candidate[0], candidate[1], groups["N1"], groups["N3"], n1,
                                    logit_scale) - \
                _logpost_gp(gp_pre[0], gp_pre[1], groups["N1"], groups["N3"], n1, logit_scale)
            if not math.isnan(log_ratio):
                if rmath.log(rng.runif(1)[0]) < log_ratio:
                    gp_pre = candidate
            par = gp_parameters(gp_pre, version)
            s["delta"] = rbase.r_min(par[0], K)
            s["nu"] = float(par[1])

        # Atom probabilities of a semiparametric or nonparametric f1.
        if code == 2:
            s["psi_sp"] = draw_atoms(u, s["psi_sp"], concentration, s, gamma, C, atoms, version)
        elif code == 3:
            s["psi_np"] = draw_atoms(u, s["psi_np"], concentration, s, gamma, C, atoms, version)

        fdr = mixture_fdr(model, u, s, C, version)["fdr"]
        if store:
            if C is None:
                raise ValueError("no cutoff could be drawn")
            kept["C"][h] = C
            for name in ("tau", "pi0", "lambda", "eta", "theta", "xi", "delta", "nu"):
                kept[name][h] = s[name]
            fdr_draws[h] = fdr
            h += 1

    first = version["fdr_from"](burn)
    if first > draws:
        raise IndexError("subscript out of bounds")
    if first == draws or U < 2:
        raise ValueError("'x' must be an array of at least two dimensions")
    return {"values": u, "frequency": frequency, "fdr": rbase.col_means(fdr_draws[first - 1:draws]),
            **kept, "acceptance": acceptance}
