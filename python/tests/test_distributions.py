# -----------------------------------------------------------------------------
# File:    python/tests/test_distributions.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Tests of distributions.py: the ZIGP probabilities against their
#          formula, the samplers against their means, and the same doubles
#          and draws as R.
# Usage:   python -m pytest python/tests/test_distributions.py
# -----------------------------------------------------------------------------
import math

import numpy as np
import pytest

from localfdr import rng
from localfdr.distributions import (check_zigp, dzigp, rdirichlet, rmultinomial, rzigp,
                                    shifted_density)


def test_check_zigp():
    check_zigp(0.4, 2.0, 0.3)
    for eta, lam, theta in ((1.2, 2.0, 0.3), (0.4, 0.0, 0.3), (0.4, 2.0, 1.5), (math.nan, 2, 0.3)):
        with pytest.raises(ValueError):
            check_zigp(eta, lam, theta)


def test_dzigp():
    x = np.arange(200)
    assert dzigp(x, 0.4, 2.0, 0.3).sum() == pytest.approx(1.0, abs=1e-12)
    # theta = 0: zero-inflated Poisson
    poisson = np.array([math.exp(-2.0) * 2.0 ** k / math.factorial(k) for k in range(6)])
    expected = 0.6 * poisson
    expected[0] += 0.4
    np.testing.assert_allclose(dzigp(np.arange(6), 0.4, 2.0, 0.0), expected, rtol=1e-14)
    assert dzigp(0, 0.25, 1.0, 0.5)[0] == pytest.approx(0.25 + 0.75 * math.exp(-1.0))
    with pytest.raises(ValueError):
        dzigp(x, 1.5, 2.0, 0.3)
    assert np.isnan(dzigp([1.0], 0.5, -1.0, 0.3, check=False)[0])


def test_rzigp():
    rng.set_seed(1)
    draws = np.array([rzigp(0.4, 2.0, 0.3) for _ in range(20000)])
    assert abs((draws == 0).mean() - (0.4 + 0.6 * math.exp(-2.0))) < 0.01
    assert abs(draws.mean() - 0.6 * 2.0 / 0.7) < 0.05   # (1 - eta) lambda / (1 - theta)
    with pytest.raises(ValueError):
        rzigp(0.4, -2.0, 0.3)


def test_rdirichlet():
    rng.set_seed(2)
    draws = np.array([rdirichlet([1.0, 2.0, 7.0]) for _ in range(5000)])
    np.testing.assert_allclose(draws.sum(axis=1), 1.0, rtol=1e-14)
    np.testing.assert_allclose(draws.mean(axis=0), [0.1, 0.2, 0.7], atol=0.01)
    assert np.isnan(rdirichlet([1.0, np.nan])).all()
    replaced = rdirichlet([1.0, np.nan], replace_missing=True)
    assert replaced.sum() == pytest.approx(1.0) and replaced[0] > 0.5


def test_rmultinomial():
    rng.set_seed(3)
    draws = np.array([rmultinomial(10, [0.2, 0.3, 0.5]) for _ in range(5000)])
    assert (draws.sum(axis=1) == 10).all()
    np.testing.assert_allclose(draws.mean(axis=0), [2.0, 3.0, 5.0], atol=0.1)
    np.testing.assert_array_equal(rmultinomial(1, [0.0, 1.0, 0.0]), [0.0, 1.0, 0.0])
    np.testing.assert_array_equal(rmultinomial(4, [0.0, 0.0, 3.0]), [0.0, 0.0, 4.0])
    with pytest.raises(ValueError):
        rmultinomial(1, [0.5, np.nan])


def test_shifted_density():
    x = np.arange(12)
    poisson = shifted_density("P", x, 2.0, 1.0, 0.5, 5)
    assert (poisson[:6] == 0).all()
    np.testing.assert_allclose(poisson[6:], [math.exp(-2.0) * 2.0 ** k / math.factorial(k)
                                             for k in range(6)], rtol=1e-14)
    generalized = shifted_density("GP", np.arange(300), 2.0, 1.5, 0.3, 5)
    assert (generalized[:6] == 0).all() and generalized.sum() == pytest.approx(1.0, abs=1e-10)
    assert (shifted_density("NP", x, 2.0, 1.0, 0.5, 5) == 0).all()


@pytest.mark.against_r
def test_same_doubles_and_draws_as_r(run_r, reader):
    r = reader(run_r(r"""
x <- 0:150
for (p in list(c(0.4, 2, 0.3), c(0, 0.01, 1), c(1, 30, 0), c(0.25, 7.5, 0.95))) put(dzigp(x, p[1], p[2], p[3]))
put(shifted_density("P", x, 2.5, 1, 0.5, 5)); put(shifted_density("GP", x, 1, 1.7, 0.4, 10))
put(suppressWarnings(shifted_density("GP", x, 1, 1.7, -0.4, 10)))
set.seed(31)
for (p in list(c(0.4, 2, 0.3), c(0.1, 0.3, 0.99), c(0, 40, 0.5), c(0.7, 4, 0))) put(replicate(400, rzigp(p[1], p[2], p[3])))
put(rzigp(0.01, 600, 0.999))
for (i in 1:50) put(rdirichlet(c(0.02, 0.5, 1, 1.5, 3, 20.5, 1e-5)))
put(suppressWarnings(rdirichlet(c(1, NA, 2)))); put(rdirichlet(c(1, NA, 2), TRUE))
q <- c(0.1, 0, 0.3, 0.05, 0, 0.55)
for (i in 1:200) put(rmultinomial(1, q)); for (i in 1:50) put(rmultinomial(37, q * 3))
put(rmultinomial(5, c(0, 0, 2))); put(runif(1))
"""))
    x = np.arange(151)
    for eta, lam, theta in ((0.4, 2, 0.3), (0, 0.01, 1), (1, 30, 0), (0.25, 7.5, 0.95)):
        r.same(f"dzigp {eta} {lam} {theta}", dzigp(x, eta, lam, theta))
    r.same("shifted Poisson", shifted_density("P", x, 2.5, 1, 0.5, 5))
    r.same("shifted GP", shifted_density("GP", x, 1, 1.7, 0.4, 10))
    r.same("shifted GP, negative dispersion", shifted_density("GP", x, 1, 1.7, -0.4, 10))
    rng.set_seed(31)
    for eta, lam, theta in ((0.4, 2, 0.3), (0.1, 0.3, 0.99), (0, 40, 0.5), (0.7, 4, 0)):
        r.same(f"rzigp {eta} {lam} {theta}", [rzigp(eta, lam, theta) for _ in range(400)])
    r.same("rzigp, probabilities on the log scale", [rzigp(0.01, 600, 0.999)])
    for _ in range(50):
        r.same("rdirichlet", rdirichlet([0.02, 0.5, 1, 1.5, 3, 20.5, 1e-5]))
    r.same("rdirichlet with a missing value", rdirichlet([1, np.nan, 2]))
    r.same("rdirichlet, missing value replaced", rdirichlet([1, np.nan, 2], True))
    q = np.array([0.1, 0, 0.3, 0.05, 0, 0.55])
    for _ in range(200):
        r.same("rmultinomial 1", rmultinomial(1, q))
    for _ in range(50):
        r.same("rmultinomial 37", rmultinomial(37, q * 3))
    r.same("rmultinomial, last category", rmultinomial(5, [0, 0, 2]))
    r.same("stream", rng.runif(1))
    r.done()
