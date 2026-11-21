# -----------------------------------------------------------------------------
# File:    python/localfdr/eb_fdr.py
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
# Usage:   from localfdr.eb_fdr import eb_rejections
#          eb_rejections(x, level=0.05)   # number of positions declared non-null
# Source:  R/eb_fdr.R. Every sum and product is rounded in the order R
#          evaluates it, so the fit is the same doubles as in R. Positions
#          are counted from zero here (`storey_rejected()`); the cutoff is a
#          count, as in R.
# -----------------------------------------------------------------------------
import math

import numpy as np
from numba import njit

from . import rbase, rmath
from .distributions import _log
from .rmath import _lgammafn

# Small constant the fit adds to probabilities and to denominators.
EB_FLOOR = math.pow(10.0, -10.0)
_EM_TOLERANCE = math.pow(10.0, -4.0)


def _vector(x):
    return np.atleast_1d(np.asarray(x, dtype=np.float64))


# ----- ZIGP null fitted by EM -----

@njit(cache=True, error_model="numpy")
def _gp(x, lam, theta, lfact):
    """Generalized Poisson probability of one count; lfact is lfactorial(x)."""
    return math.exp(_log(lam) + (x - 1) * _log(lam + theta * x) - lfact - lam - theta * x)


@njit(cache=True, error_model="numpy")
def _gp_pmf(x, lam, theta):
    out = np.empty(len(x))
    for i in range(len(x)):
        out[i] = _gp(x[i], lam, theta, _lgammafn(x[i] + 1))
    return out


def gp_pmf(x, lam, theta):
    """Generalized Poisson probabilities.

    Parameters: x (counts), lam, theta (rate and dispersion).
    Returns: array of probabilities.
    """
    return _gp_pmf(_vector(x), float(lam), float(theta))


def zigp_pmf_floored(x, eta, lam, theta):
    """ZIGP probabilities, kept away from zero.

    Parameters: x (counts), eta, lam, theta (ZIGP parameters).
    Returns: array of probabilities.
    """
    x = _vector(x)
    return eta * (1 * (x == 0)) + (1 - eta) * (gp_pmf(x, lam, theta) + EB_FLOOR) + EB_FLOOR


@njit(cache=True, error_model="numpy")
def _gp_sum(C, lam, theta, lfact):
    """sum(gp_pmf(0:C, lam, theta))."""
    total = 0.0
    for v in range(C + 1):
        total += _gp(float(v), lam, theta, lfact[v])
    return total


def zigp_cdf(eta, lam, theta, C):
    """Probability that a ZIGP count is at most C.

    Parameters: eta, lam, theta (ZIGP parameters), C (cutoff).
    Returns: the probability (float).
    """
    return float(eta + (1 - eta) * rbase.r_sum(gp_pmf(np.arange(int(C) + 1), lam, theta)))


@njit(cache=True, error_model="numpy")
def _em_tail(K, eta, lam, theta, C, n, lfact):
    """(n, total, total_rate, n_rate, total_dispersion, n_dispersion) of
    em_tail(); lfact[v] is lfactorial(v) for v = 0, ..., K + 1."""
    inside = _gp_sum(C, lam, theta, lfact)
    below = eta + (1 - eta) * inside
    total = total_rate = n_rate = total_dispersion = n_dispersion = 0.0
    count = K - C if C < K else 2  # (C + 1):K runs downwards when C = K
    for step in range(count):
        v = float(C + 1 + step) if C < K else float(C + 1 - step)
        pmf = _gp(v, lam, theta, lfact[int(v)])
        pmf_times_v = math.exp(_log(v * lam) + (v - 1) * _log(lam + theta * v) - lfact[int(v)]
                               - lam - theta * v)
        rate = lam / (lam + theta * v)
        dispersion = (theta * v) / (lam + theta * v)
        total += pmf_times_v
        total_rate += pmf_times_v * rate
        n_rate += pmf * rate
        total_dispersion += pmf_times_v * dispersion
        n_dispersion += pmf * dispersion
    return (n * ((1 - eta) * (1 - inside)) / below,
            n * ((1 - eta) * total) / below,
            n * ((1 - eta) * total_rate) / below,
            n * ((1 - eta) * n_rate) / below,
            n * ((1 - eta) * total_dispersion) / below,
            n * ((1 - eta) * n_dispersion) / below)


def em_tail(K, eta, lam, theta, C, n):
    """Expected contributions of the counts above the cutoff in the E-step.

    The fit sees only the n counts up to C; the counts above C are replaced
    by their expectation under the current ZIGP, summed up to the largest
    count K.

    Parameters: K (largest observed count), eta, lam, theta (current ZIGP
    parameters), C (cutoff), n (number of counts up to C).
    Returns: dict of the expected number of counts above C (`n`), their sum
    (`total`), and the sums weighted by lambda / (lambda + theta x) and by
    theta x / (lambda + theta x) (`total_rate`, `n_rate`,
    `total_dispersion`, `n_dispersion`).
    """
    K, C = int(K), int(C)
    lfact = rmath.lfactorial(np.arange(K + 2))
    values = _em_tail(K, float(eta), float(lam), float(theta), C, float(n), lfact)
    names = ("n", "total", "total_rate", "n_rate", "total_dispersion", "n_dispersion")
    return dict(zip(names, (float(v) for v in values)))


@njit(cache=True, error_model="numpy")
def _floor_at(lower, x):
    """max(lower, x), NaN when x is missing."""
    if math.isnan(x):
        return x
    return lower if x < lower else x


@njit(cache=True, error_model="numpy")
def _zigp_em(count, K, eta_start, lam_start, theta_start, lfact, floor, epsilon, estimates):
    """Fill estimates (3 x K) with the EM estimates for every cutoff; returns
    False when an estimate is missing (R stops there)."""
    for C in range(1, K + 1):
        n = 0.0
        for i in range(C + 1):
            n += count[i]
        eta = eta_start
        lam = lam_start
        theta = theta_start
        while True:
            eta0 = _floor_at(0.0, eta)
            lam0 = _floor_at(floor, lam)
            theta0 = _floor_at(floor, theta)
            tail_n, total, total_rate, n_rate, total_dispersion, n_dispersion = \
                _em_tail(K, eta0, lam0, theta0, C, n, lfact)
            u1 = u2 = weighted = u6 = v1 = 0.0
            for i in range(C + 1):
                value = float(i)
                lagged = float(i - 1) if i >= 1 else 0.0
                gp = _gp(value, lam0, theta0, lfact[i])
                mixture = eta0 * (1.0 if i == 0 else 0.0) + (1 - eta0) * (gp + floor) + floor
                from_gp = (1 - eta0) * (gp + floor) / mixture
                if i == 0:
                    u1 = count[0] * (eta0 * 1.0 / mixture)
                u2 += count[i] * from_gp
                weighted += value * count[i] * from_gp
                u6 += lagged * count[i] * (lam0 / (lam0 + theta0 * value))
                v1 += lagged * count[i] * ((theta0 * value) / (lam0 + theta0 * value))
            eta = u1 / (u1 + u2 + tail_n + floor)
            A = weighted + total
            B = u2 + tail_n + u6 + total_rate - n_rate
            lam = (B + floor) / (u2 + tail_n + floor)
            D = v1 + total_dispersion - n_dispersion
            theta = (D + floor) / (A + floor)
            moved = abs(eta0 - eta) > epsilon or abs(lam0 - lam) > epsilon or \
                abs(theta0 - theta) > epsilon
            if not moved:
                if math.isnan(eta) or math.isnan(lam) or math.isnan(theta) or \
                        math.isnan(eta0) or math.isnan(lam0) or math.isnan(theta0):
                    return False
                break
        estimates[0, C - 1] = eta
        estimates[1, C - 1] = lam
        estimates[2, C - 1] = theta
    return True


def zigp_em(x):
    """EM estimates of the ZIGP null for every cutoff C = 1, ..., max(x).

    For each C the ZIGP is fitted to the counts up to C, treating the counts
    above C as missing. Starting values are moment estimates from the lower
    three quarters of the data; the iterations stop when no parameter moves
    by more than 1e-4.

    Parameters: x (counts; for an array of integers the mean of the lower
    three quarters is the sum divided by the length, as R's mean() of an
    integer vector, for an array of floats it has R's second pass).
    Returns: matrix with rows eta, lambda, theta and one column per cutoff.
    """
    whole = np.issubdtype(np.asarray(x).dtype, np.integer)  # an integer vector in R
    x = _vector(x)
    K = int(x.max())
    lower = x[x <= rbase.quantile(x, 0.75)]
    ex = rbase.r_mean_integer(lower) if whole else rbase.r_mean(lower)
    vx = rbase.r_var(lower)
    zero_share = rbase.r_mean_integer(x == 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        er = float(np.float64(zero_share) / (1 - zero_share))
        inside = np.float64(vx) / np.float64(ex) - er * ex + EB_FLOOR
    phi = math.sqrt(inside) if inside >= 0 else float("nan")
    if ex != 0 and len(lower) < 2:
        raise ValueError("missing value where TRUE/FALSE needed")
    if ex == 0 or phi == 0 or math.isnan(phi):
        theta_start = 1.0
        lam_start = 0.0
    else:
        theta_start = min(max(0.0, 1 - 1 / phi), 1.0)
        lam_start = max(0.0, ex / ((1 - zero_share) * phi))
    eta_start = max(0.0, zero_share)

    count = rbase.count_frequencies_from_zero(x, K)
    estimates = np.zeros((3, K))
    lfact = rmath.lfactorial(np.arange(K + 2))
    if not _zigp_em(count, K, eta_start, lam_start, theta_start, lfact, EB_FLOOR, _EM_TOLERANCE,
                    estimates):
        raise ValueError("missing value where TRUE/FALSE needed")
    return estimates


# ----- Cutoff and null proportion -----

def zigp_negloglik(count, eta, lam, theta, C):
    """Negative log likelihood of the ZIGP over the counts up to C.

    Parameters: count (number of observations at 0, 1, ..., max(x)), eta,
    lam, theta (ZIGP parameters), C (cutoff).
    Returns: the negative log likelihood (float).
    """
    count = _vector(count)
    C = int(C)
    v = np.arange(1, C + 1, dtype=np.float64)
    at_zero = count[0] * (-rmath.log(eta + (1 - eta) * rmath.exp(-lam)))
    positive = count[1:C + 1] * (-rmath.log((1 - eta) * gp_pmf(v, lam, theta)))
    return float(at_zero + rbase.r_sum(positive))


def eb_fit(x):
    """Fit the ZIGP null, the null proportion and the cutoff by empirical
    Bayes.

    For each candidate cutoff C the ZIGP is fitted by EM and pi0 is the share
    of counts up to C divided by the fitted probability of that range. The
    cutoff is the first local minimum (or else the minimum) of the negative
    log likelihood of the counts up to C plus that of the counts above C
    under their observed frequencies.

    Parameters: x (counts).
    Returns: dict: `estimates` (matrix from `zigp_em()`), `pi0` (per
    cutoff), `criterion` (per cutoff), `cutoff`, and `eta`, `lambda`,
    `theta`, `pi0_hat` at the chosen cutoff.
    """
    estimates = zigp_em(x)
    x = _vector(x)
    K = int(x.max())
    M = len(x)
    count = rbase.count_frequencies_from_zero(x, K)
    pi0 = np.zeros(K)
    criterion = np.zeros(K)
    log_M = rmath.log(M)
    for C in range(1, K + 1):
        eta, lam, theta = estimates[:, C - 1]
        with np.errstate(divide="ignore", invalid="ignore"):
            share = np.float64(rbase.r_sum(count[:C + 1])) / (M * zigp_cdf(eta, lam, theta, C))
        pi0[C - 1] = rbase.r_min(1, share)
        if C < K:
            tail = count[C + 1:K + 1]
            above = log_M * float(np.sum(x > C)) - rbase.r_sum(tail * rmath.log(tail + EB_FLOOR))
        else:
            above = log_M * float(np.sum(x > C)) - count[K - 1] * rmath.log(count[K - 1] + EB_FLOOR)
        criterion[C - 1] = above + zigp_negloglik(count, eta, lam, theta, C)
    with np.errstate(invalid="ignore"):
        turns = np.flatnonzero(np.diff(np.sign(np.diff(criterion))) > 0) + 3
    if np.isnan(criterion).all():
        raise ValueError("no cutoff: the criterion is missing for every cutoff")
    cutoff = int(min([int(np.nanargmin(criterion)) + 1] + [int(turn) for turn in turns]))
    return {"estimates": estimates, "pi0": pi0, "criterion": criterion, "cutoff": cutoff,
            "eta": float(estimates[0, cutoff - 1]), "lambda": float(estimates[1, cutoff - 1]),
            "theta": float(estimates[2, cutoff - 1]), "pi0_hat": float(pi0[cutoff - 1])}


def eb_local_fdr(fit, x, inverse_frequency):
    """Empirical Bayes local fdr of sorted counts.

    fdr(x) = pi0 f0(x) / f(x), with f(x) the observed relative frequency of x.

    Parameters: fit (output of `eb_fit()`), x (counts at which to evaluate),
    inverse_frequency (number of observations divided by the number of them
    equal to each element of `x`).
    Returns: array of local fdr values.
    """
    x = _vector(x)
    return fit["pi0_hat"] * (fit["eta"] * (1 * (x == 0)) +
                             (1 - fit["eta"]) * (gp_pmf(x, fit["lambda"], fit["theta"]) + EB_FLOOR) +
                             EB_FLOOR) * _vector(inverse_frequency)


def _below_level(fdr, level):
    """sum(fdr[!is.na(fdr)] < level)."""
    return float(np.sum(fdr[~np.isnan(fdr)] < level))


def eb_rejections(x, level=0.05):
    """Number of positions declared non-null in a data set (one-stage
    procedure).

    Parameters: x (counts, one per position), level (fdr level).
    Returns: number of positions with local fdr below `level` (float).
    """
    fit = eb_fit(x)
    x = _vector(x)
    values, frequency = np.unique(x, return_counts=True)
    ordered = np.sort(x)
    fdr = eb_local_fdr(fit, ordered, len(x) / frequency[np.searchsorted(values, ordered)])
    return _below_level(fdr, level)


# ----- The three procedures in a simulated data set -----

def second_stage_threshold(fit, n):
    """Second-stage threshold: counts from here on are rejected outright.

    Parameters: fit (output of `eb_fit()`), n (number of counts).
    Returns: the threshold (float).
    """
    lam, theta = fit["lambda"], fit["theta"]
    with np.errstate(divide="ignore", invalid="ignore"):
        first = float(np.float64(lam) / (rmath.exp(theta - 1) - theta))
        second = float(np.float64(rmath.log(n)) / (theta - 1 - rmath.log(theta)))
    return float(np.ceil(rbase.r_max(fit["cutoff"] + 1, rbase.r_max(first, second))))


def storey_rejected(data, fit, level):
    """Rejections of Storey's procedure with ZIGP tail p-values.

    The p-value of a count is its ZIGP upper tail probability. The sorted
    p-values are compared with (number of counts at least as large as each
    distinct value) * level / (n pi0).

    Parameters: data (output of `simulated_data()`), fit (output of
    `eb_fit()`), level (FDR level).
    Returns: positions (among the distinct counts, from zero) that are
    rejected.
    """
    x = _vector(data["x"])
    n = len(x)
    largest = int(x.max())
    pmf = gp_pmf(np.arange(largest), fit["lambda"], fit["theta"])
    tail = np.empty(largest)
    inside = 0.0
    for i in range(largest):
        inside += pmf[i]
        tail[i] = rbase.r_max(0, (1 - fit["eta"]) * (1 - inside))
    pvalue = np.concatenate(([1.0], tail))[_vector(data["values"]).astype(np.int64)]
    frequency = _vector(data["null_frequency"]) + _vector(data["alt_frequency"])
    if len(frequency) < 2:
        raise ValueError("a data set needs two distinct counts")
    by_pvalue = rbase.r_order(pvalue)
    at_least = np.empty(len(frequency))
    at_least[0] = n
    running = 0.0
    for i in range(1, len(frequency)):
        running += frequency[i - 1]
        at_least[i] = n - running
    with np.errstate(divide="ignore", invalid="ignore"):
        bound = at_least * (np.float64(level) / (n * fit["pi0_hat"]))
    if np.isnan(pvalue).any() or np.isnan(bound).any():
        raise ValueError("missing value where TRUE/FALSE needed")
    below = pvalue[by_pvalue] <= bound
    if below.sum() == 0:
        return np.zeros(0, dtype=np.int64)
    last = int(np.flatnonzero(below).max()) + 1
    if last == 1:
        return np.zeros(0, dtype=np.int64)  # a single rejection is not kept
    return by_pvalue[:last]


def eb_simulation_counts(data, C_true, level):
    """Rejection counts of the two-stage, one-stage and Storey procedures in
    one simulated data set.

    Two-stage: null counts below the second-stage threshold are rejected by
    their local fdr, every count from the threshold on is rejected.
    One-stage: counts are rejected by their local fdr alone.

    Parameters: data (output of `simulated_data()`), C_true (cutoff the data
    were simulated with), level (fdr level).
    Returns: dict (one row) with TS_R, TS_V, TS_TP, OS_R, OS_V, OS_TP,
    Storey_R, Storey_V, Storey_TP (rejections, false and true rejections).
    """
    null = _vector(data["null"])
    alt = _vector(data["alt"])
    values = _vector(data["values"])
    null_frequency = _vector(data["null_frequency"])
    alt_frequency = _vector(data["alt_frequency"])
    n = len(data["x"])
    fit = eb_fit(data["x"])
    threshold = second_stage_threshold(fit, n)

    null_sorted = np.sort(null)
    alt_sorted = np.sort(alt)

    def weight(x):
        return n / (null_frequency + alt_frequency)[np.searchsorted(values, x)]

    alt_weight = weight(alt_sorted)

    def below_level(x, weights=None):
        return _below_level(eb_local_fdr(fit, x, weight(x) if weights is None else weights), level)

    # Null counts rejected by their local fdr.
    null_below_threshold = below_level(null_sorted[null_sorted < threshold])
    null_up_to_cutoff = below_level(null_sorted[null_sorted <= fit["cutoff"]])
    null_above_cutoff = below_level(null_sorted[null_sorted > fit["cutoff"]])
    # Non-null counts rejected by their local fdr; the weights are taken from
    # the smallest non-null counts in both ranges.
    n_alt_below = int(np.sum(alt < threshold))
    n_alt_from = int(np.sum(alt > threshold - 1))
    alt_below_threshold = below_level(alt_sorted[alt_sorted < threshold], alt_weight[:n_alt_below])
    alt_from_threshold = below_level(alt_sorted[alt_sorted > (threshold - 1)], alt_weight[:n_alt_from])

    TS_V = max(0.0, null_below_threshold - null_up_to_cutoff) + float(np.sum(null >= threshold))
    TS_TP = alt_below_threshold + float(np.sum(alt >= threshold))
    OS_V = max(0.0, null_above_cutoff)
    OS_TP = alt_below_threshold + alt_from_threshold

    storey = storey_rejected(data, fit, level)
    Storey_V = rbase.r_sum(null_frequency[storey])
    Storey_TP = rbase.r_sum(alt_frequency[storey])

    return {"TS_R": TS_TP + TS_V, "TS_V": TS_V, "TS_TP": TS_TP,
            "OS_R": OS_TP + OS_V, "OS_V": OS_V, "OS_TP": OS_TP,
            "Storey_R": Storey_V + Storey_TP, "Storey_V": Storey_V, "Storey_TP": Storey_TP}
