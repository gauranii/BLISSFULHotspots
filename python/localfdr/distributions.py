# -----------------------------------------------------------------------------
# File:    python/localfdr/distributions.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Tables 1 to 3, Figures S1 to S3.
# Purpose: The zero-inflated generalized Poisson (ZIGP) null distribution, the
#          shifted Poisson and generalized Poisson alternatives, and the
#          Dirichlet and multinomial samplers of the Gibbs sampler.
# Inputs:  None.
# Outputs: None.
# Usage:   from localfdr.distributions import dzigp, rzigp, rdirichlet
# Source:  R/distributions.R. Every sum and product is rounded in the order
#          R evaluates it, and the draws use the random numbers of rng.py,
#          so the results are the same doubles as in R.
# -----------------------------------------------------------------------------
import math

import numpy as np
from numba import njit

from . import rng
from .rmath import _lgammafn
from .rng import _rbinom, _rgamma, _unif


def check_zigp(eta, lam, theta):
    """Stop when a ZIGP parameter is outside its range (or missing).

    Parameters: eta, lam, theta. Returns: None.
    """
    if not (eta >= 0 and eta <= 1):
        raise ValueError("eta >= 0 && eta <= 1 is not TRUE")
    if not lam > 0:
        raise ValueError("lambda > 0 is not TRUE")
    if not (theta >= 0 and theta <= 1):
        raise ValueError("theta >= 0 && theta <= 1 is not TRUE")


@njit(cache=True, error_model="numpy")
def _log(x):
    return math.log(x) if x >= 0.0 else np.nan


@njit(cache=True, error_model="numpy")
def _dzigp(x, eta, lam, theta):
    density = np.zeros(len(x))
    log_lam = _log(lam)
    for i in range(len(x)):
        gp = math.exp(log_lam + (x[i] - 1) * _log(lam + theta * x[i]) - _lgammafn(x[i] + 1)
                      - (lam + theta * x[i]))
        if x[i] == 0:
            density[i] = eta + (1 - eta) * gp
        elif x[i] != 0:  # false for a missing count, whose density stays 0
            density[i] = (1 - eta) * gp
    return density


def dzigp(x, eta, lam, theta, check=True):
    """Probability mass function of the zero-inflated generalized Poisson.

    P(X = 0) = eta + (1 - eta) exp(-lambda) and, for x > 0,
    P(X = x) = (1 - eta) lambda (lambda + theta x)^(x - 1) exp(-lambda - theta x) / x!.
    theta = 0 gives the zero-inflated Poisson.

    Parameters: x (counts, non-negative), eta (zero-inflation probability,
    in [0, 1]), lam (rate, positive), theta (dispersion, in [0, 1]), check
    (if True stop when a parameter is outside its range).
    Returns: array of probabilities, as long as `x`.
    """
    if check:
        check_zigp(eta, lam, theta)
    return _dzigp(np.atleast_1d(np.asarray(x, dtype=np.float64)), float(eta), float(lam),
                  float(theta))


@njit(cache=True, error_model="numpy")
def _rzigp(state, eta, lam, theta):
    """One ZIGP draw by inversion; -1 when a probability is missing (R
    stops there)."""
    unif = _unif(state)
    cdf = eta + (1 - eta) * math.exp(-lam)
    if math.isnan(cdf):
        return -1.0
    recursive = False
    log_p_previous = 0.0
    p = 0.0
    j = 1
    while cdf < unif:
        if not recursive:
            p = math.exp(_log(lam) + (j - 1) * _log(lam + theta * j) - _lgammafn(j + 1.0)
                         - (lam + theta * j) + _log(1 - eta))
        if math.isnan(p):
            return -1.0
        if p == np.inf:
            recursive = True
            log_p_previous = _log(lam) + (j - 2) * _log(lam + theta * (j - 1)) - _lgammafn(float(j)) \
                - (lam + theta * (j - 1)) + _log(1 - eta)
        if recursive:
            log_p = log_p_previous + _log(lam + theta * (j - 1)) - _log(float(j)) - theta \
                + (j - 1) * _log(lam + theta * j) - (j - 1) * _log(lam + theta * (j - 1))
            log_p_previous = log_p
            p = math.exp(log_p)
        cdf = cdf + p
        if math.isnan(cdf):
            return -1.0
        j += 1
    return float(j - 1)


def rzigp(eta, lam, theta, check=True):
    """One draw from the zero-inflated generalized Poisson, by inversion.

    Uses one uniform number. The probabilities are summed from zero upwards
    until they pass the uniform number; if a probability overflows it is
    continued on the log scale from the previous one.

    Parameters: eta, lam, theta (as in `dzigp()`), check (if True stop when
    a parameter is outside its range).
    Returns: one count (float).
    """
    if check:
        check_zigp(eta, lam, theta)
    draw = _rzigp(rng._state, float(eta), float(lam), float(theta))
    if draw < 0:
        raise ValueError("missing value where TRUE/FALSE needed")
    return float(draw)


@njit(cache=True, error_model="numpy")
def _rdirichlet(state, beta, replace_missing):
    draw = np.empty(len(beta))
    for i in range(len(beta)):
        shape = beta[i]
        if replace_missing and math.isnan(shape):
            shape = 0.00001
        draw[i] = _rgamma(state, shape, 1.0)
    total = 0.0
    for i in range(len(beta)):  # summed in order, in double precision
        total = total + draw[i]
    return draw / total


def rdirichlet(beta, replace_missing=False):
    """One draw from a Dirichlet distribution, from independent gamma draws.

    Parameters: beta (concentration parameters), replace_missing (if True a
    missing concentration is set to 1e-5).
    Returns: probabilities that sum to one, as long as `beta`.
    """
    return _rdirichlet(rng._state, np.atleast_1d(np.asarray(beta, dtype=np.float64)),
                       bool(replace_missing))


@njit(cache=True, error_model="numpy")
def _rmultinomial(state, size, prob, draw):
    """Fill `draw`; returns False where R stops (a missing probability or
    number of trials)."""
    K = len(prob)
    for i in range(K):
        draw[i] = 0.0
    remaining = 0.0
    for i in range(K):
        remaining += prob[i]
    for i in range(K):
        if math.isnan(prob[i]):
            return False
        if prob[i] != 0:
            p = prob[i] / remaining
            if math.isnan(p):
                return False
            draw[i] = _rbinom(state, size, p) if p < 1 else size
            size = size - draw[i]
        if math.isnan(size):
            return False
        if size <= 0:
            return True
        remaining = remaining - prob[i]
    if K >= 2:
        draw[K - 2] = size
    return True


def rmultinomial(size, prob):
    """One draw from a multinomial distribution, as a chain of binomial draws.

    Parameters: size (number of trials), prob (probabilities of the
    categories; need not sum to one).
    Returns: counts per category (float array).
    """
    prob = np.atleast_1d(np.asarray(prob, dtype=np.float64))
    draw = np.empty(len(prob))
    if not _rmultinomial(rng._state, float(size), prob, draw):
        raise ValueError("missing value where TRUE/FALSE needed")
    return draw


@njit(cache=True, error_model="numpy")
def _shifted_density(poisson, x, xi, delta, nu, C, f1):
    """Fill f1 with the shifted Poisson (poisson True) or generalized
    Poisson density."""
    for i in range(len(x)):
        f1[i] = 0.0
        if x[i] > C:
            shift = x[i] - C - 1
            if poisson:
                f1[i] = math.exp(shift * _log(xi) - _lgammafn(shift + 1) - xi)
            else:
                f1[i] = math.exp(_log(delta) - _lgammafn(shift + 1)
                                 + (shift - 1) * _log(delta + nu * shift) - (delta + nu * shift))


def shifted_density(model, x, xi, delta, nu, C):
    """Parametric alternative density f1: a Poisson or generalized Poisson
    shifted to start right above the cutoff C.

    Parameters: model ("P": Poisson with mean `xi`; "GP": generalized
    Poisson with rate `delta` and dispersion `nu`), x (counts), xi, delta,
    nu (parameters of f1), C (cutoff: f1 is zero for counts up to C).
    Returns: array of f1(x).
    """
    x = np.atleast_1d(np.asarray(x, dtype=np.float64))
    f1 = np.zeros(len(x))
    if model in ("P", "GP"):
        _shifted_density(model == "P", x, float(xi), float(delta), float(nu), float(C), f1)
    return f1
