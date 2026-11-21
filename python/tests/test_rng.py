# -----------------------------------------------------------------------------
# File:    python/tests/test_rng.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Tests of rng.py: known draws of R, the range and the mean of every
#          sampler, and the same streams as R for several seeds and
#          parameters (every branch of the gamma, beta, binomial and Poisson
#          algorithms).
# Usage:   python -m pytest python/tests/test_rng.py
# -----------------------------------------------------------------------------
import numpy as np
import pytest

from localfdr import rng

SEEDS = (1234, 1, 889192, 1440, 2147483647)


def test_set_seed_gives_the_draws_of_r():
    rng.set_seed(1)
    np.testing.assert_allclose(rng.runif(3), [0.2655087, 0.3721239, 0.5728534], atol=5e-8)
    rng.set_seed(1)
    np.testing.assert_allclose(rng.rnorm(3), [-0.6264538, 0.1836433, -0.8356286], atol=5e-8)
    rng.set_seed(1)
    np.testing.assert_array_equal(rng.sample_int(5), [1, 4, 3, 5, 2])
    rng.set_seed(42)
    first = rng.runif(4)
    rng.set_seed(42)
    np.testing.assert_array_equal(rng.runif(4), first)


def test_rnorm_and_rexp_moments():
    rng.set_seed(3)
    z = rng.rnorm(40000, 2.0, 0.5)
    assert abs(z.mean() - 2.0) < 0.01 and abs(z.std() - 0.5) < 0.01
    e = rng.rexp(40000)
    assert e.min() > 0 and abs(e.mean() - 1.0) < 0.02


def test_rgamma():
    rng.set_seed(4)
    for shape, rate in ((0.5, 1.0), (2.5, 2.0), (20.0, 0.5)):
        draws = rng.rgamma(40000, shape, rate)
        assert draws.min() > 0
        assert abs(draws.mean() - shape / rate) < 0.03 * shape / rate
    out = rng.rgamma(4, [0.0, np.nan, -1.0, 2.0], 1.0)
    assert out[0] == 0.0 and np.isnan(out[1]) and np.isnan(out[2]) and out[3] > 0
    assert len(rng.rgamma(0, 1.0)) == 0


def test_rbeta():
    rng.set_seed(5)
    for a, b in ((0.5, 0.5), (1.0, 3.0), (6.0, 2.0)):
        draws = rng.rbeta(40000, a, b)
        assert draws.min() >= 0 and draws.max() <= 1
        assert abs(draws.mean() - a / (a + b)) < 0.01
    assert np.isnan(rng.rbeta(1, np.nan, 1.0)[0]) and np.isnan(rng.rbeta(1, -1.0, 1.0)[0])
    assert rng.rbeta(1, 0.0, 2.0)[0] == 0.0 and rng.rbeta(1, 2.0, 0.0)[0] == 1.0


def test_rbinom():
    rng.set_seed(6)
    for size, prob in ((1, 0.3), (20, 0.9), (250, 0.2), (1000, 0.5)):
        draws = rng.rbinom(40000, size, prob)
        assert draws.min() >= 0 and draws.max() <= size and (draws == np.rint(draws)).all()
        assert abs(draws.mean() - size * prob) < 0.02 * size * prob + 0.01
    rng.set_seed(6)
    next_uniform = rng.runif(1)[0]
    rng.set_seed(6)
    out = rng.rbinom(4, 7, [0.0, 1.0, np.nan, 1.5])  # none of these uses the stream
    assert out[0] == 0 and out[1] == 7 and np.isnan(out[2]) and np.isnan(out[3])
    assert rng.runif(1)[0] == next_uniform


def test_rpois_and_rgeom():
    rng.set_seed(7)
    for mean in (0.4, 5.0, 12.0, 300.0):
        draws = rng.rpois(40000, mean)
        assert draws.min() >= 0 and (draws == np.rint(draws)).all()
        assert abs(draws.mean() - mean) < 0.02 * mean + 0.01
    assert rng.rpois(1, 0.0)[0] == 0 and np.isnan(rng.rpois(1, -1.0)[0])
    geometric = rng.rgeom(40000, 0.25)
    assert abs(geometric.mean() - 3.0) < 0.08
    assert np.isnan(rng.rgeom(1, 0.0)[0])


def test_sample_int_and_sample():
    rng.set_seed(8)
    assert sorted(rng.sample_int(20)) == list(range(1, 21))
    assert len(set(rng.sample_int(100, 30))) == 30
    with_replacement = rng.sample_int(3, 3000, replace=True)
    assert set(with_replacement) == {1, 2, 3}
    weighted = rng.sample_int(3, 30000, replace=True, prob=[1, 2, 7])
    np.testing.assert_allclose(np.bincount(weighted)[1:] / 30000, [0.1, 0.2, 0.7], atol=0.01)
    walker = rng.sample_int(300, 60000, replace=True, prob=np.arange(1, 301))
    assert abs(walker.mean() - 200.33) < 1.5
    assert sorted(rng.sample_int(4, 4, prob=[1, 2, 3, 4])) == [1, 2, 3, 4]
    assert set(rng.sample_int(3, 2, prob=[1, 0, 1])) == {1, 3}
    letters = np.array(list("abcde"))
    assert sorted(rng.sample(letters)) == list("abcde")
    with pytest.raises(ValueError):
        rng.sample_int(3, 4)
    with pytest.raises(ValueError):
        rng.sample_int(3, 2, prob=[1, 2])


def python_streams(parameters=None):
    """The draws of R_STREAMS from rng.py, as (label, values) pairs."""
    out = []
    for seed in SEEDS:
        rng.set_seed(seed)
        out.append(("runif", rng.runif(5)))
        out.append(("rnorm", rng.rnorm(7)))
        out.append(("rexp", rng.rexp(50)))
        out.append(("rgamma", rng.rgamma(300, [0.02, 0.5, 0.999, 1, 1.5, 2.7, 3.686, 5, 13.022, 14,
                                               80.5, 1e-3], [1, 0.5, 1 + 1e-10])))
        out.append(("rbeta", rng.rbeta(300, [1, 0.3, 2, 30, 151.5, 1.01], [1, 4, 0.7, 2, 45, 1.02])))
        out.append(("rbinom", rng.rbinom(400, [1, 1, 5, 250, 50, 1000, 3000],
                                         [0.3, 0.97, 0.2, 0.5, 0.811, 0.031, 0.5, 0.999])))
        out.append(("rpois", rng.rpois(300, [0.3, 2, 5, 9.99, 10, 12.5, 40, 1000])))
        out.append(("rgeom", rng.rgeom(100, [0.25, 0.01, 0.9])))
        out.append(("sample.int", rng.sample_int(10)))
        out.append(("sample.int replace", rng.sample_int(7, 20, replace=True)))
        out.append(("sample.int prob replace",
                    rng.sample_int(6, 30, replace=True, prob=[0.1, 0.3, 0.05, 0.3, 0.2, 0.05])))
        out.append(("sample.int prob", rng.sample_int(6, 4, prob=[1, 3, 0.5, 3, 2, 0.5])))
        out.append(("sample.int Walker",
                    rng.sample_int(400, 50, replace=True, prob=np.tile([1, 2, 3, 4], 100))))
        out.append(("runif after", rng.runif(2)))
    return out


R_STREAMS = r"""
for (seed in c(1234, 1, 889192, 1440, 2147483647)) {
  set.seed(seed)
  put(runif(5)); put(rnorm(7)); put(rexp(50))
  put(rgamma(300, c(0.02, 0.5, 0.999, 1, 1.5, 2.7, 3.686, 5, 13.022, 14, 80.5, 1e-3), c(1, 0.5, 1 + 1e-10)))
  put(rbeta(300, c(1, 0.3, 2, 30, 151.5, 1.01), c(1, 4, 0.7, 2, 45, 1.02)))
  put(rbinom(400, c(1, 1, 5, 250, 50, 1000, 3000), c(0.3, 0.97, 0.2, 0.5, 0.811, 0.031, 0.5, 0.999)))
  put(rpois(300, c(0.3, 2, 5, 9.99, 10, 12.5, 40, 1000)))
  put(rgeom(100, c(0.25, 0.01, 0.9)))
  put(sample.int(10)); put(sample.int(7, 20, replace = TRUE))
  put(sample.int(6, 30, replace = TRUE, prob = c(0.1, 0.3, 0.05, 0.3, 0.2, 0.05)))
  put(sample.int(6, 4, prob = c(1, 3, 0.5, 3, 2, 0.5)))
  put(sample.int(400, 50, replace = TRUE, prob = rep(c(1, 2, 3, 4), 100)))
  put(runif(2))
}
"""


@pytest.mark.against_r
def test_streams_of_every_sampler_equal_r(run_r, reader):
    r = reader(run_r(R_STREAMS))
    for label, values in python_streams():
        r.same(label, values)
    r.done()


@pytest.mark.against_r
def test_long_streams_over_random_parameters_equal_r(run_r, reader):
    r = reader(run_r(r"""
set.seed(77)
shape <- c(exp(runif(400, -6, 6)), 1:40 + 0.5); rate <- exp(runif(37, -2, 2))
a <- exp(runif(301, -3, 6)); b <- exp(runif(299, -3, 6))
size <- floor(exp(runif(211, 0, 9))); prob <- runif(199); mu <- exp(runif(307, -3, 8))
put(shape); put(rate); put(a); put(b); put(size); put(prob); put(mu)
set.seed(4321)
put(rgamma(60000, shape, rate)); put(rbeta(60000, a, b)); put(rbinom(60000, size, prob))
put(rpois(60000, mu)); put(rexp(20000)); put(rnorm(20000)); put(runif(3))
"""))
    shape, rate, a, b = r.take(440), r.take(37), r.take(301), r.take(299)
    size, prob, mu = r.take(211), r.take(199), r.take(307)
    rng.set_seed(4321)
    r.same("rgamma", rng.rgamma(60000, shape, rate))
    r.same("rbeta", rng.rbeta(60000, a, b))
    r.same("rbinom", rng.rbinom(60000, size, prob))
    r.same("rpois", rng.rpois(60000, mu))
    r.same("rexp", rng.rexp(20000))
    r.same("rnorm", rng.rnorm(20000))
    r.same("runif", rng.runif(3))
    r.done()
