# -----------------------------------------------------------------------------
# File:    python/localfdr/rng.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Tables 1 to 3, Figure S2.
# Purpose: The random numbers of R (>= 3.6.0, default kinds): set.seed() with
#          its initial scrambling, the Mersenne-Twister, unif_rand(), and the
#          samplers built on it, ported operation by operation: runif(),
#          rnorm() by inversion, exponential draws (exp_rand), rgamma()
#          (Ahrens and Dieter GD for shape >= 1, GS below 1), rbeta() (Cheng
#          BB and BC), rbinom() (inversion and BTPE), rpois() (Ahrens and
#          Dieter PD), rgeom() and sample.int() (rejection sampling; with
#          probabilities: inversion, Walker's alias method, or without
#          replacement). A seed gives the same draws as in R, so the same
#          simulated data and the same Gibbs chains. One stream per process,
#          as in R.
# Inputs:  None.
# Outputs: None.
# Usage:   from localfdr import rng
#          rng.set_seed(1234); rng.runif(3); rng.rgamma(2, [0.5, 3.0], 2.0)
# Source:  R 4.6.1: src/main/RNG.c, src/main/random.c (do_sample),
#          src/main/sort.c (revsort), src/nmath/snorm.c, qnorm.c, sexp.c,
#          rgamma.c, rbeta.c, rbinom.c, rpois.c, rgeom.c, runif.c;
#          src/library/stats/src/random.c (recycling of the parameters,
#          missing values).
# Notes:   R's samplers keep the constants of the last parameters between
#          calls; these are functions of the parameters alone, so they are
#          computed at every call here. muladd() marks the products that R's
#          C code adds with one rounding on arm64 (see rmath.py).
# -----------------------------------------------------------------------------
import math

import numpy as np
from numba import njit

from .rmath import DBL_MAX, M_1_SQRT_2PI, M_LN2, _pow_di, muladd

N = 624  # words of the Mersenne-Twister state
M = 397

# The stream: N state words, then the position of the next word.
_state = np.zeros(N + 1, dtype=np.uint32)
_state[N] = N + 1



@njit(cache=True)
def _seed(state, seed):
    """Fill the state from a seed as `set.seed()` does.

    Parameters: state (uint32 array of length N + 1), seed (integer).
    Returns: None; `state` is overwritten.
    """
    s = np.uint32(seed)
    for _ in range(50):  # initial scrambling
        s = np.uint32(69069) * s + np.uint32(1)
    for j in range(N + 1):
        s = np.uint32(69069) * s + np.uint32(1)
        if j > 0:
            state[j - 1] = s
    state[N] = N  # the first word is replaced by the position


@njit(cache=True)
def _unif(state):
    """One uniform number in (0, 1): `unif_rand()`.

    Parameters: state (the stream). Returns: float.
    """
    mti = state[N]
    if mti >= N:
        for kk in range(N):
            y = (state[kk] & np.uint32(0x80000000)) | (state[(kk + 1) % N] & np.uint32(0x7FFFFFFF))
            v = state[(kk + M) % N] ^ (y >> np.uint32(1))
            if y & np.uint32(1):
                v ^= np.uint32(0x9908B0DF)
            state[kk] = v
        mti = 0
    y = state[mti]
    state[N] = mti + 1
    y ^= y >> np.uint32(11)
    y ^= (y << np.uint32(7)) & np.uint32(0x9D2C5680)
    y ^= (y << np.uint32(15)) & np.uint32(0xEFC60000)
    y ^= y >> np.uint32(18)
    x = float(y) * 2.3283064365386963e-10
    if x <= 0.0:
        return 0.5 * 2.328306437080797e-10
    if 1.0 - x <= 0.0:
        return 1.0 - 0.5 * 2.328306437080797e-10
    return x


@njit(cache=True)
def _unif_index(state, dn):
    """A uniform index in 0, ..., dn - 1 by rejection: `R_unif_index()`.

    Parameters: state (the stream), dn (number of values). Returns: int.
    """
    if dn <= 0:
        return 0
    bits = int(np.ceil(np.log2(float(dn))))
    while True:
        v = 0
        n = 0
        while n <= bits:
            v = 65536 * v + int(np.floor(_unif(state) * 65536))
            n += 16
        if bits < 64:
            v &= (1 << bits) - 1
        if v < dn:
            return v


@njit(cache=True)
def _sample_int(state, n, size):
    """`sample.int(n, size)` without replacement, values 1 to n.

    Parameters: state (the stream), n, size. Returns: int64 array.
    """
    pool = np.arange(n)
    out = np.empty(size, dtype=np.int64)
    for i in range(size):
        j = _unif_index(state, n)
        out[i] = pool[j] + 1
        n -= 1
        pool[j] = pool[n]
    return out


@njit(cache=True)
def _horner(r, coefficients):
    """Polynomial in r, highest power first, evaluated as qnorm() does.

    Parameters: r, coefficients (array).
    Returns: float.
    """
    value = coefficients[0]
    for c in coefficients[1:]:
        value = muladd(value, r, c)
    return value


_CENTRE_NUM = np.array([2509.0809287301226727, 33430.575583588128105, 67265.770927008700853,
                        45921.953931549871457, 13731.693765509461125, 1971.5909503065514427,
                        133.14166789178437745, 3.387132872796366608])
_CENTRE_DEN = np.array([5226.495278852854561, 28729.085735721942674, 39307.89580009271061,
                        21213.794301586595867, 5394.1960214247511077, 687.1870074920579083,
                        42.313330701600911252, 1.0])
_MIDDLE_NUM = np.array([7.7454501427834140764e-4, .0227238449892691845833, .24178072517745061177,
                        1.27045825245236838258, 3.64784832476320460504, 5.7694972214606914055,
                        4.6303378461565452959, 1.42343711074968357734])
_MIDDLE_DEN = np.array([1.05075007164441684324e-9, 5.475938084995344946e-4,
                        .0151986665636164571966, .14810397642748007459, .68976733498510000455,
                        1.6763848301838038494, 2.05319162663775882187, 1.0])
_TAIL_NUM = np.array([2.01033439929228813265e-7, 2.71155556874348757815e-5,
                      .0012426609473880784386, .026532189526576123093, .29656057182850489123,
                      1.7848265399172913358, 5.4637849111641143699, 6.6579046435011037772])
_TAIL_DEN = np.array([2.04426310338993978564e-15, 1.4215117583164458887e-7,
                      1.8463183175100546818e-5, 7.868691311456132591e-4,
                      .0148753612908506148525, .13692988092273580531, .59983220655588793769,
                      1.0])


@njit(cache=True)
def _qnorm(p):
    """Standard normal quantile of a probability in (0, 1): `qnorm()`.

    Parameters: p. Returns: float.
    """
    q = p - 0.5
    if abs(q) <= 0.425:
        r = muladd(-q, q, 0.180625)
        return q * _horner(r, _CENTRE_NUM) / _horner(r, _CENTRE_DEN)
    r = np.sqrt(-np.log(1.0 - p if q > 0 else p))
    if r <= 5.0:
        r += -1.6
        value = _horner(r, _MIDDLE_NUM) / _horner(r, _MIDDLE_DEN)
    else:
        r += -5.0
        value = _horner(r, _TAIL_NUM) / _horner(r, _TAIL_DEN)
    return -value if q < 0.0 else value


# ----- Normal and exponential draws -----

@njit(cache=True)
def _norm_rand(state):
    """One standard normal number by inversion: `norm_rand()`.

    Two uniform numbers per value. Parameters: state. Returns: float.
    """
    big = 134217728.0  # 2^27
    u = _unif(state)
    u = int(big * u) + _unif(state)
    return _qnorm(u / big)


_EXP_Q = np.array([0.6931471805599453, 0.9333736875190459, 0.9888777961838675,
                   0.9984959252914960, 0.9998292811061389, 0.9999833164100727,
                   0.9999985691438767, 0.9999998906925558, 0.9999999924734159,
                   0.9999999995283275, 0.9999999999728814, 0.9999999999985598,
                   0.9999999999999289, 0.9999999999999968, 0.9999999999999999,
                   1.0000000000000000])


@njit(cache=True)
def _exp_rand(state):
    """One standard exponential number: `exp_rand()` (Ahrens and Dieter 1972).

    Parameters: state. Returns: float.
    """
    q = _EXP_Q
    a = 0.0
    u = _unif(state)
    while u <= 0.0 or u >= 1.0:
        u = _unif(state)
    while True:
        u += u
        if u > 1.0:
            break
        a += q[0]
    u -= 1.0
    if u <= q[0]:
        return a + u
    i = 0
    ustar = _unif(state)
    umin = ustar
    while True:
        ustar = _unif(state)
        if umin > ustar:
            umin = ustar
        i += 1
        if not u > q[i]:
            break
    return muladd(umin, q[0], a)


# ----- Gamma -----

@njit(cache=True, error_model="numpy")
def _gamma_quotient(t, s, s2, q0):
    """The quotient q of steps 6 and 10 of algorithm GD."""
    v = t / (s + s)
    if abs(v) <= 0.25:
        poly = muladd(0.1233795, v, -0.1367177)
        poly = muladd(poly, v, 0.1423657)
        poly = muladd(poly, v, -0.1662921)
        poly = muladd(poly, v, 0.2000062)
        poly = muladd(poly, v, -0.250003)
        poly = muladd(poly, v, 0.3333333)
        return muladd(0.5 * t * t * poly, v, q0)
    q = muladd(-s, t, q0)
    q = muladd(0.25 * t, t, q)
    return muladd(s2 + s2, math.log(1.0 + v), q)


@njit(cache=True, error_model="numpy")
def _rgamma(state, a, scale):
    """One gamma number with shape a and this scale: `rgamma()` of the C
    library. NaN for a parameter outside its range.

    Parameters: state, a, scale. Returns: float.
    """
    if math.isnan(a) or math.isnan(scale):
        return np.nan
    if a <= 0.0 or scale <= 0.0:
        if scale == 0.0 or a == 0.0:
            return 0.0
        return np.nan
    if math.isinf(a) or math.isinf(scale):
        return np.inf

    if a < 1.0:  # algorithm GS
        e = muladd(0.36787944117144232159, a, 1.0)
        while True:
            p = e * _unif(state)
            if p >= 1.0:
                x = -math.log((e - p) / a)
                if _exp_rand(state) >= (1.0 - a) * math.log(x):
                    break
            else:
                x = math.exp(math.log(p) / a)
                if _exp_rand(state) >= x:
                    break
        return scale * x

    # algorithm GD
    s2 = a - 0.5
    s = math.sqrt(s2)
    d = muladd(-s, 12.0, 5.656854)
    t = _norm_rand(state)
    x = muladd(0.5, t, s)
    ret_val = x * x
    if t >= 0.0:
        return scale * ret_val
    u = _unif(state)
    if d * u <= t * t * t:
        return scale * ret_val

    r = 1.0 / a
    q0 = muladd(2.424e-4, r, 2.4511e-4)
    q0 = muladd(q0, r, -7.388e-5)
    q0 = muladd(q0, r, 0.00144121)
    q0 = muladd(q0, r, 0.00801191)
    q0 = muladd(q0, r, 0.02083148)
    q0 = muladd(q0, r, 0.04166669)
    q0 = q0 * r
    if a <= 3.686:
        b = muladd(0.178, s2, 0.463 + s)
        si = 1.235
        c = muladd(0.16, s, 0.195 / s - 0.079)
    elif a <= 13.022:
        b = muladd(0.0076, s2, 1.654)
        si = 1.68 / s + 0.275
        c = 0.062 / s + 0.024
    else:
        b = 1.77
        si = 0.75
        c = 0.1515 / s

    if x > 0.0:
        q = _gamma_quotient(t, s, s2, q0)
        if math.log(1.0 - u) <= q:
            return scale * ret_val

    while True:
        e = _exp_rand(state)
        u = _unif(state)
        u = u + u - 1.0
        if u < 0.0:
            t = muladd(-si, e, b)
        else:
            t = muladd(si, e, b)
        if t >= -0.71874483771719:
            q = _gamma_quotient(t, s, s2, q0)
            if q > 0.0:
                w = math.expm1(q)
                if c * abs(u) <= w * math.exp(muladd(-(0.5 * t), t, e)):
                    break
    x = muladd(0.5, t, s)
    return scale * x * x


# ----- Beta -----

_EXPMAX = 1024 * M_LN2  # log(DBL_MAX)


@njit(cache=True, error_model="numpy")
def _beta_w(beta, u1, factor):
    """v and w of Cheng's algorithms from the first uniform number."""
    v = beta * math.log(u1 / (1.0 - u1))
    if v <= _EXPMAX:
        w = factor * math.exp(v)
        if math.isinf(w):
            w = DBL_MAX
    else:
        w = DBL_MAX
    return v, w


@njit(cache=True, error_model="numpy")
def _rbeta(state, aa, bb):
    """One beta number: `rbeta()` of the C library (Cheng 1978, BB and BC).

    Parameters: state, aa, bb (shapes). Returns: float, NaN for a shape that
    is missing or negative.
    """
    if math.isnan(aa) or math.isnan(bb) or aa < 0.0 or bb < 0.0:
        return np.nan
    if math.isinf(aa) and math.isinf(bb):
        return 0.5
    if aa == 0.0 and bb == 0.0:
        return 0.0 if _unif(state) < 0.5 else 1.0
    if math.isinf(aa) or bb == 0.0:
        return 1.0
    if math.isinf(bb) or aa == 0.0:
        return 0.0

    a = min(aa, bb)
    b = max(aa, bb)
    alpha = a + b
    if a <= 1.0:  # algorithm BC
        beta = 1.0 / a
        delta = 1.0 + b - a
        k1 = delta * muladd(0.0416667, a, 0.0138889) / muladd(b, beta, -0.777778)
        k2 = muladd(0.5 + 0.25 / delta, a, 0.25)
        while True:
            u1 = _unif(state)
            u2 = _unif(state)
            if u1 < 0.5:
                y = u1 * u2
                z = u1 * y
                if muladd(0.25, u2, z) - y >= k1:
                    continue
            else:
                z = u1 * u1 * u2
                if z <= 0.25:
                    v, w = _beta_w(beta, u1, b)
                    break
                if z >= k2:
                    continue
            v, w = _beta_w(beta, u1, b)
            if muladd(alpha, math.log(alpha / (a + w)) + v, -1.3862944) >= math.log(z):
                break
        return a / (a + w) if aa == a else w / (a + w)

    # algorithm BB
    beta = math.sqrt((alpha - 2.0) / muladd(2.0 * a, b, -alpha))
    gamma = a + 1.0 / beta
    while True:
        u1 = _unif(state)
        u2 = _unif(state)
        v, w = _beta_w(beta, u1, a)
        z = u1 * u1 * u2
        r = muladd(gamma, v, -1.3862944)
        s = a + r - w
        if s + 2.609438 >= 5.0 * z:
            break
        t = math.log(z)
        if s > t:
            break
        if not muladd(alpha, math.log(alpha / (b + w)), r) < t:
            break
    return b / (b + w) if aa != a else w / (b + w)


# ----- Binomial -----

@njit(cache=True, error_model="numpy")
def _stirling_tail(x2, x1):
    return (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0


@njit(cache=True, error_model="numpy")
def _rbinom(state, nin, pp):
    """One binomial number: `rbinom()` of the C library (inversion when
    n p < 30, else BTPE of Kachitvichyanukul and Schmeiser 1988).

    Parameters: state, nin (number of trials), pp (probability).
    Returns: float, NaN for a parameter that is missing or out of range.
    """
    if math.isnan(nin) or math.isinf(nin):
        return np.nan
    r = np.rint(nin)
    if r != nin:
        return np.nan
    if math.isnan(pp) or math.isinf(pp) or r < 0 or pp < 0.0 or pp > 1.0:
        return np.nan
    if r == 0 or pp == 0.0:
        return 0.0
    if pp == 1.0:
        return r
    if r >= 2147483647:
        return np.nan  # R draws by qbinom() here; not ported
    n = int(r)

    p = min(pp, 1.0 - pp)
    q = 1.0 - p
    np_ = n * p
    r = p / q
    g = r * (n + 1)
    ix = 0

    if np_ < 30.0:  # inversion
        qn = _pow_di(q, n)
        while True:
            ix = 0
            f = qn
            u = _unif(state)
            found = False
            while True:
                if u < f:
                    found = True
                    break
                if ix > 110:
                    break
                u -= f
                ix += 1
                f *= (g / ix - r)
            if found:
                break
    else:  # BTPE
        ffm = np_ + p
        m = int(ffm)
        fm = float(m)
        npq = np_ * q
        p1 = int(muladd(2.195, math.sqrt(npq), -(4.6 * q))) + 0.5
        xm = fm + 0.5
        xl = xm - p1
        xr = xm + p1
        c = 0.134 + 20.5 / (15.3 + fm)
        al = (ffm - xl) / muladd(-xl, p, ffm)
        xll = al * muladd(0.5, al, 1.0)
        al = (xr - ffm) / (xr * q)
        xlr = al * muladd(0.5, al, 1.0)
        p2 = p1 * (1.0 + c + c)
        p3 = p2 + c / xll
        p4 = p3 + c / xlr
        while True:
            u = _unif(state) * p4
            v = _unif(state)
            if u <= p1:  # triangular region
                ix = int(muladd(-p1, v, xm) + u)
                break
            if u <= p2:  # parallelogram region
                x = xl + (u - p1) / c
                v = muladd(v, c, 1.0) - abs(xm - x) / p1
                if v > 1.0 or v <= 0.0:
                    continue
                ix = int(x)
            elif u > p3:  # right tail
                ix = int(xr - math.log(v) / xlr)
                if ix > n:
                    continue
                v = v * (u - p3) * xlr
            else:  # left tail
                ix = int(xl + math.log(v) / xll)
                if ix < 0:
                    continue
                v = v * (u - p2) * xll
            k = abs(ix - m)
            if k <= 20 or k >= npq / 2 - 1:  # explicit evaluation
                f = 1.0
                if m < ix:
                    for i in range(m + 1, ix + 1):
                        f *= (g / i - r)
                elif m != ix:
                    for i in range(ix + 1, m + 1):
                        f /= (g / i - r)
                if v <= f:
                    break
            else:  # squeezing with bounds on log(f(x))
                amaxp = (k / npq) * (muladd(float(k), k / 3.0 + 0.625, 0.1666666666666) / npq + 0.5)
                ynorm = -k * k / (2.0 * npq)
                alv = math.log(v)
                if alv < ynorm - amaxp:
                    break
                if alv <= ynorm + amaxp:
                    x1 = float(ix + 1)
                    f1 = fm + 1.0
                    z = n + 1 - fm
                    w = n - ix + 1.0
                    z2 = z * z
                    x2 = x1 * x1
                    f2 = f1 * f1
                    w2 = w * w
                    bound = muladd(xm, math.log(f1 / x1), (n - m + 0.5) * math.log(z / w))
                    bound = muladd(float(ix - m), math.log(w * p / (x1 * q)), bound)
                    bound = bound + _stirling_tail(f2, f1) + _stirling_tail(z2, z) + \
                        _stirling_tail(x2, x1) + _stirling_tail(w2, w)
                    if alv <= bound:
                        break
    if pp > 0.5:
        ix = n - ix
    return float(ix)


# ----- Poisson and geometric -----

_FACT = np.array([1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0])


@njit(cache=True, error_model="numpy")
def _rpois(state, mu):
    """One Poisson number: `rpois()` of the C library (Ahrens and Dieter
    1982: table look-up below a mean of 10, normal deviates from 10 on).

    Parameters: state, mu (mean). Returns: float, NaN for a mean that is
    missing, infinite or negative.
    """
    if math.isnan(mu) or math.isinf(mu) or mu < 0:
        return np.nan
    if mu <= 0.0:
        return 0.0

    if mu < 10.0:  # inversion over a table of cumulative probabilities
        m = max(1, int(mu))
        table_end = 0  # the table pp[] is filled up to this index
        pp = np.empty(36)
        p0 = math.exp(-mu)
        p = p0
        q = p0
        while True:
            u = _unif(state)
            if u <= p0:
                return 0.0
            if table_end != 0:
                start = 1 if u <= 0.458 else min(table_end, m)
                for k in range(start, table_end + 1):
                    if u <= pp[k]:
                        return float(k)
                if table_end == 35:
                    continue
            table_end += 1
            for k in range(table_end, 36):
                p *= mu / k
                q += p
                pp[k] = q
                if u <= q:
                    return float(k)
            table_end = 35

    s = math.sqrt(mu)
    d = 6.0 * mu * mu
    big_l = math.floor(mu - 1.1484)
    g = muladd(s, _norm_rand(state), mu)
    pois = -1.0
    u = 0.0
    fk = 0.0
    difmuk = 0.0
    if g >= 0.0:
        pois = math.floor(g)
        if pois >= big_l:
            return pois
        fk = pois
        difmuk = mu - fk
        u = _unif(state)
        if d * u >= difmuk * difmuk * difmuk:
            return pois

    omega = M_1_SQRT_2PI / s
    b1 = 0.0416666666666666667 / mu
    b2 = 0.3 * b1 * b1
    c3 = 0.1428571428571428571 * b1 * b2
    c2 = muladd(-15.0, c3, b2)
    c1 = muladd(45.0, c3, muladd(-6.0, b2, b1))
    c0 = muladd(-15.0, c3, muladd(3.0, b2, 1.0 - b1))
    c = 0.1069 / mu

    quotient_step = g >= 0.0  # first pass with the normal deviate (step Q)
    E = 0.0
    while True:
        if not quotient_step:
            E = _exp_rand(state)
            u = muladd(2.0, _unif(state), -1.0)
            t = 1.8 + (abs(E) if u >= 0 else -abs(E))
            if not t > -0.6744:
                continue
            pois = math.floor(muladd(s, t, mu))
            fk = pois
            difmuk = mu - fk
        # 'subroutine' F: px, py, fx, fy
        if pois < 10:
            px = -mu
            py = math.pow(mu, pois) / _FACT[int(pois)]
        else:
            delta = 0.0833333333333333333 / fk
            delta = delta * muladd(-(4.8 * delta), delta, 1.0)
            v = difmuk / fk
            if abs(v) <= 0.25:
                poly = muladd(0.1250060, v, -0.1384794)
                poly = muladd(poly, v, 0.1421878)
                poly = muladd(poly, v, -0.1661269)
                poly = muladd(poly, v, 0.2000118)
                poly = muladd(poly, v, -0.2500068)
                poly = muladd(poly, v, 0.3333333)
                poly = muladd(poly, v, -0.5)
                px = muladd(fk * v * v, poly, -delta)
            else:
                px = muladd(fk, math.log(1.0 + v), -difmuk) - delta
            py = M_1_SQRT_2PI / math.sqrt(fk)
        x = (0.5 - difmuk) / s
        x *= x
        fx = -0.5 * x
        fy = omega * muladd(muladd(muladd(c3, x, c2), x, c1), x, c0)
        if not quotient_step:
            if c * abs(u) <= muladd(py, math.exp(px + E), -(fy * math.exp(fx + E))):
                break
        else:
            quotient_step = False
            if muladd(-u, fy, fy) <= py * math.exp(px - fx):
                break
    return pois


@njit(cache=True, error_model="numpy")
def _rgeom(state, p):
    """One geometric number: `rgeom()` of the C library.

    Parameters: state, p (success probability). Returns: float.
    """
    if math.isnan(p) or math.isinf(p) or p <= 0 or p > 1:
        return np.nan
    return _rpois(state, _exp_rand(state) * ((1 - p) / p))


# ----- Vectors of draws, parameters recycled as in R -----

@njit(cache=True)
def _runif(state, n):
    out = np.empty(n)
    for i in range(n):
        out[i] = _unif(state)
    return out


@njit(cache=True)
def _rnorm(state, n, mean, sd):
    out = np.empty(n)
    for i in range(n):
        out[i] = muladd(sd, _norm_rand(state), mean)
    return out


@njit(cache=True)
def _rexp_std(state, n):
    out = np.empty(n)
    for i in range(n):
        out[i] = _exp_rand(state)
    return out


@njit(cache=True)
def _rgamma_n(state, n, shape, scale):
    out = np.empty(n)
    for i in range(n):
        out[i] = _rgamma(state, shape[i % len(shape)], scale[i % len(scale)])
    return out


@njit(cache=True)
def _rbeta_n(state, n, shape1, shape2):
    out = np.empty(n)
    for i in range(n):
        out[i] = _rbeta(state, shape1[i % len(shape1)], shape2[i % len(shape2)])
    return out


@njit(cache=True)
def _rbinom_n(state, n, size, prob):
    out = np.empty(n)
    for i in range(n):
        out[i] = _rbinom(state, size[i % len(size)], prob[i % len(prob)])
    return out


@njit(cache=True)
def _rpois_n(state, n, mean):
    out = np.empty(n)
    for i in range(n):
        out[i] = _rpois(state, mean[i % len(mean)])
    return out


@njit(cache=True)
def _rgeom_n(state, n, prob):
    out = np.empty(n)
    for i in range(n):
        out[i] = _rgeom(state, prob[i % len(prob)])
    return out


def _vector(x):
    return np.atleast_1d(np.asarray(x, dtype=np.float64))


def set_seed(seed):
    """`set.seed(seed)`: start the stream of this process at a seed.

    Parameters: seed (integer). Returns: None.
    """
    _seed(_state, int(seed) & 0xFFFFFFFF)


def runif(n):
    """`runif(n)`: n uniform numbers in (0, 1). Returns: float array."""
    return _runif(_state, n)


def rnorm(n, mean=0.0, sd=1.0):
    """`rnorm(n, mean, sd)`: n normal numbers. Returns: float array."""
    return _rnorm(_state, n, float(mean), float(sd))


def rexp(n):
    """`rexp(n)`: n standard exponential numbers. Returns: float array."""
    return _rexp_std(_state, n)


def rgamma(n, shape, rate=1.0):
    """`rgamma(n, shape, rate)`: n gamma numbers; shape and rate are recycled.

    Parameters: n, shape, rate (floats or arrays).
    Returns: float array, NaN where a parameter is missing or out of range.
    """
    return _rgamma_n(_state, n, _vector(shape), 1.0 / _vector(rate))


def rbeta(n, shape1, shape2):
    """`rbeta(n, shape1, shape2)`: n beta numbers. Returns: float array."""
    return _rbeta_n(_state, n, _vector(shape1), _vector(shape2))


def rbinom(n, size, prob):
    """`rbinom(n, size, prob)`: n binomial numbers; size and prob are recycled.

    Returns: float array, NaN (R: NA) where a parameter is missing or out
    of range.
    """
    return _rbinom_n(_state, n, _vector(size), _vector(prob))


def rpois(n, mean):
    """`rpois(n, mean)`: n Poisson numbers. Returns: float array."""
    return _rpois_n(_state, n, _vector(mean))


def rgeom(n, prob):
    """`rgeom(n, prob)`: n geometric numbers. Returns: float array."""
    return _rgeom_n(_state, n, _vector(prob))


# ----- sample.int() -----

@njit(cache=True)
def _sample_replace(state, n, size):
    out = np.empty(size, dtype=np.int64)
    for i in range(size):
        out[i] = _unif_index(state, n) + 1
    return out


@njit(cache=True)
def _revsort(a, ib):
    """Sort a into descending order by heapsort, ib alongside: revsort()."""
    n = len(a)
    if n <= 1:
        return
    left = (n >> 1) + 1
    ir = n
    while True:
        if left > 1:
            left -= 1
            ra = a[left - 1]
            ii = ib[left - 1]
        else:
            ra = a[ir - 1]
            ii = ib[ir - 1]
            a[ir - 1] = a[0]
            ib[ir - 1] = ib[0]
            ir -= 1
            if ir == 1:
                a[0] = ra
                ib[0] = ii
                return
        i = left
        j = left << 1
        while j <= ir:
            if j < ir and a[j - 1] > a[j]:
                j += 1
            if ra > a[j - 1]:
                a[i - 1] = a[j - 1]
                ib[i - 1] = ib[j - 1]
                i = j
                j += i
            else:
                j = ir + 1
        a[i - 1] = ra
        ib[i - 1] = ii


@njit(cache=True)
def _sample_prob(state, p, size, replace):
    """sample.int(n, size, replace, prob) for normalized probabilities p
    (overwritten): ProbSampleReplace(), walker_ProbSampleReplace() and
    ProbSampleNoReplace()."""
    n = len(p)
    out = np.empty(size, dtype=np.int64)
    perm = np.arange(1, n + 1)
    if replace:
        large = 0
        for i in range(n):
            if n * p[i] > 0.1:
                large += 1
        if large > 200:  # Walker's alias method
            hl = np.empty(n, dtype=np.int64)
            q = np.empty(n)
            alias = np.zeros(n, dtype=np.int64)
            low = 0   # hl[0:low]: entries with q < 1
            high = n  # hl[high:n]: entries with q >= 1
            for i in range(n):
                q[i] = p[i] * n
                if q[i] < 1.0:
                    hl[low] = i
                    low += 1
                else:
                    high -= 1
                    hl[high] = i
            if low > 0 and high < n:
                for k in range(n - 1):
                    i = hl[k]
                    j = hl[high]
                    alias[i] = j
                    q[j] += q[i] - 1.0
                    if q[j] < 1.0:
                        high += 1
                    if high >= n:
                        break
            for i in range(n):
                q[i] += i
            for i in range(size):
                k = _unif_index(state, n)
                u = k + _unif(state)
                out[i] = k + 1 if u < q[k] else alias[k] + 1
            return out
        _revsort(p, perm)
        for i in range(1, n):
            p[i] += p[i - 1]
        for i in range(size):
            u = _unif(state)
            j = 0
            while j < n - 1:
                if u <= p[j]:
                    break
                j += 1
            out[i] = perm[j]
        return out
    _revsort(p, perm)
    total = 1.0
    n1 = n - 1
    for i in range(size):
        target = total * _unif(state)
        mass = 0.0
        j = 0
        while j < n1:
            mass += p[j]
            if target <= mass:
                break
            j += 1
        out[i] = perm[j]
        total -= p[j]
        for k in range(j, n1):
            p[k] = p[k + 1]
            perm[k] = perm[k + 1]
        n1 -= 1
    return out


def sample_int(n, size=None, replace=False, prob=None):
    """`sample.int(n, size, replace, prob)`: `size` of 1, ..., n.

    Parameters: n, size (n when None), replace (with replacement), prob
    (weights of 1, ..., n, or None for equal weights).
    Returns: int64 array, values 1 to n.
    """
    size = n if size is None else size
    if not replace and size > n:
        raise ValueError("cannot take a sample larger than the population when 'replace = FALSE'")
    if prob is None:
        if n > 1e7 and not replace and size <= n / 2:
            raise NotImplementedError("sampling by hashing (n > 1e7) is not ported")
        if replace or size < 2:
            return _sample_replace(_state, n, size)
        return _sample_int(_state, n, size)
    p = np.array(prob, dtype=np.float64)
    if len(p) != n:
        raise ValueError("incorrect number of probabilities")
    if not np.all(np.isfinite(p)) or np.any(p < 0):
        raise ValueError("NA in probability vector" if np.any(np.isnan(p))
                         else "negative probability")
    positive = int(np.sum(p > 0))
    if positive == 0 or (not replace and size > positive):
        raise ValueError("too few positive probabilities")
    total = 0.0
    for value in p:
        total += value
    p = p / total
    return _sample_prob(_state, p, size, bool(replace))


def sample(x, size=None, replace=False, prob=None):
    """`sample(x, size, replace, prob)` for a vector of length above one.

    Parameters: x (array), size, replace, prob (as `sample_int()`).
    Returns: array of the sampled elements.
    """
    x = np.asarray(x)
    return x[sample_int(len(x), size, replace, prob) - 1]
