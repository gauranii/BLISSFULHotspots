# -----------------------------------------------------------------------------
# File:    python/tests/test_rmath.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Tests of rmath.py: known values of every function, and the same
#          doubles as R over grids of arguments.
# Usage:   python -m pytest python/tests/test_rmath.py
# -----------------------------------------------------------------------------
import math

import numpy as np
import pytest
from numba import njit

from localfdr import rmath


def test_muladd_adds_a_product_with_one_rounding():
    @njit
    def fused(a, b, c):
        return rmath.muladd(a, b, c)

    # 0.1 * 10 is 1 after rounding; the exact product is 1 + 5.55e-17.
    expected = 5.551115123125783e-17 if rmath.FUSED else 0.0
    assert fused(0.1, 10.0, -1.0) == expected
    assert fused(2.0, 3.0, 4.0) == 10.0


def test_lgamma_and_lfactorial_known_values():
    assert rmath.lgamma(1.0) == 0.0
    assert rmath.lgamma(5.0) == pytest.approx(math.log(24), rel=1e-15)
    assert rmath.lgamma(0.5) == pytest.approx(0.5 * math.log(math.pi), rel=1e-15)
    assert rmath.lgamma(100.0) == pytest.approx(math.lgamma(100.0), rel=1e-15)
    assert rmath.lfactorial(3) == pytest.approx(math.log(6), rel=1e-15)
    np.testing.assert_allclose(rmath.lgamma(np.array([2.0, 3.0, 11.0])),
                               [0.0, math.log(2), math.log(3628800)], rtol=1e-15)
    np.testing.assert_allclose(rmath.lfactorial(np.array([0, 1, 4])), [0.0, 0.0, math.log(24)],
                               atol=1e-15)


def test_log_and_exp_answer_every_double():
    assert rmath.log(0.0) == -math.inf
    assert math.isnan(rmath.log(-1.0))
    assert rmath.log(math.e) == 1.0
    assert rmath.exp(1000.0) == math.inf
    assert rmath.exp(0.0) == 1.0
    np.testing.assert_array_equal(rmath.log(np.array([1.0, 0.0])), [0.0, -np.inf])
    np.testing.assert_array_equal(rmath.exp(np.array([0.0, -np.inf])), [1.0, 0.0])


def test_logistic_and_uniform():
    assert rmath.plogis(0.0) == 0.5
    assert rmath.plogis(math.inf) == 1.0 and rmath.plogis(-math.inf) == 0.0
    assert rmath.plogis(1.0, 0, 0.5) == pytest.approx(1 / (1 + math.exp(-2)), rel=1e-15)
    assert rmath.dlogis(0.0) == 0.25
    assert rmath.dlogis(0.0, 0, 0.5) == 0.5
    assert math.isnan(rmath.plogis(math.nan)) and math.isnan(rmath.dlogis(1.0, 0, -1))
    assert rmath.dunif_log(0.3) == 0.0
    assert rmath.dunif_log(0.2, 0, 0.5) == pytest.approx(math.log(2))
    assert rmath.dunif_log(0.7, 0, 0.5) == -math.inf
    assert math.isnan(rmath.dunif_log(math.nan))


def test_pow_di_and_fround():
    assert rmath.pow_di(2.0, 10) == 1024.0
    assert rmath.pow_di(2.0, -2) == 0.25
    assert rmath.pow_di(3.0, 0) == 1.0
    assert rmath.fround(2.675, 2) == 2.67   # 2.675 is stored a little below
    assert rmath.fround(0.5) == 0.0 and rmath.fround(1.5) == 2.0 and rmath.fround(-2.5) == -2.0
    assert rmath.fround(0.15, 1) == 0.1 and rmath.fround(0.25, 1) == 0.2
    assert rmath.fround(1234.5678, -2) == 1200.0
    assert rmath.fround(-0.04235, 4) == -0.0424 or rmath.fround(-0.04235, 4) == -0.0423
    assert math.isnan(rmath.fround(math.nan, 2)) and rmath.fround(math.inf, 2) == math.inf


@pytest.mark.against_r
def test_same_doubles_as_r(run_r, reader):
    r = reader(run_r(r"""
set.seed(21); x <- rnorm(20000, 0, 3); put(x)
put(exp(x)); put(log(abs(x))); put(plogis(x)); put(plogis(x, 0, 0.5)); put(dlogis(x))
put(dlogis(x, 0, 0.5)); put(lgamma(abs(x) * 20)); put(lfactorial(0:200))
put(dunif(plogis(x), 0, 0.5, log = TRUE))
y <- c(x, round(x, 3) + 0.0005, 0.5, 1.5, 2.5, -0.5, 2.675, 1.005, 0.285, 1e15 + 0.5, 1e-300)
put(y); for (d in c(0, 1, 2, 4, -1)) put(round(y, d))
put(10^(-10)); put(10^(-4))
"""))
    x = r.take(20000)
    r.same("exp", rmath.exp(x))
    r.same("log", rmath.log(np.abs(x)))
    r.same("plogis", [rmath.plogis(v) for v in x])
    r.same("plogis scale", [rmath.plogis(v, 0, 0.5) for v in x])
    r.same("dlogis", [rmath.dlogis(v) for v in x])
    r.same("dlogis scale", [rmath.dlogis(v, 0, 0.5) for v in x])
    r.same("lgamma", rmath.lgamma(np.abs(x) * 20))
    r.same("lfactorial", rmath.lfactorial(np.arange(201)))
    r.same("dunif", [rmath.dunif_log(rmath.plogis(v), 0, 0.5) for v in x])
    y = r.take(40009)
    for digits in (0, 1, 2, 4, -1):
        r.same(f"round {digits}", [rmath.fround(v, digits) for v in y])
    r.same("powers", [math.pow(10.0, -10.0), math.pow(10.0, -4.0)])
    r.done()
