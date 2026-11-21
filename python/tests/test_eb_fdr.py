# -----------------------------------------------------------------------------
# File:    python/tests/test_eb_fdr.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Tests of eb_fdr.py: the generalized Poisson probabilities, the EM
#          fit, the cutoff and the three procedures on small inputs, and the
#          same doubles as R for the fit and the rejection counts.
# Usage:   python -m pytest python/tests/test_eb_fdr.py
# -----------------------------------------------------------------------------
import math

import numpy as np
import pytest

from inputs import COUNTS, SETTING
from localfdr import eb_fdr as e
from localfdr import rng
from localfdr.distributions import dzigp, rzigp
from localfdr.simulate import simulated_data


def test_floor_and_vector():
    assert e.EB_FLOOR == 1e-10 and e._vector(2).tolist() == [2.0]
    assert e._below_level(np.array([0.01, np.nan, 0.2, 0.04]), 0.05) == 2.0


def test_probabilities():
    poisson = [math.exp(-2.0) * 2.0 ** k / math.factorial(k) for k in range(8)]
    np.testing.assert_allclose(e.gp_pmf(np.arange(8), 2.0, 0.0), poisson, rtol=1e-14)
    assert e.gp_pmf(np.arange(400), 2.0, 0.4).sum() == pytest.approx(1.0, abs=1e-12)
    floored = e.zigp_pmf_floored([0, 1, 3], 0.4, 2.0, 0.3)
    np.testing.assert_allclose(floored, dzigp([0, 1, 3], 0.4, 2.0, 0.3) + (1 - 0.4) * 1e-10 + 1e-10,
                               rtol=1e-14)
    assert e.zigp_cdf(0.4, 2.0, 0.3, 3) == pytest.approx(dzigp(np.arange(4), 0.4, 2.0, 0.3).sum())
    assert e.zigp_cdf(0.4, 2.0, 0.3, 300) == pytest.approx(1.0)


def test_em_tail():
    tail = e.em_tail(60, 0.4, 2.0, 0.3, 4, 100)
    below = e.zigp_cdf(0.4, 2.0, 0.3, 4)
    v = np.arange(5, 61)
    pmf = e.gp_pmf(v, 2.0, 0.3)
    assert tail["n"] == pytest.approx(100 * (1 - below) / below, rel=1e-9)
    assert tail["total"] == pytest.approx(100 * 0.6 * (v * pmf).sum() / below, rel=1e-12)
    assert tail["n_rate"] == pytest.approx(100 * 0.6 * (pmf * 2 / (2 + 0.3 * v)).sum() / below)
    assert tail["total_rate"] + tail["total_dispersion"] == pytest.approx(tail["total"])
    assert tail["n_rate"] + tail["n_dispersion"] == pytest.approx(100 * 0.6 * pmf.sum() / below)
    assert set(e.em_tail(5, 0.4, 2.0, 0.3, 5, 10)) == set(tail)   # cutoff at the largest count


def test_zigp_em_fits_the_share_of_zeros():
    rng.set_seed(5)
    x = np.array([rzigp(0.4, 2.0, 0.3) for _ in range(4000)])
    estimates = e.zigp_em(x)
    assert estimates.shape == (3, int(x.max()))
    eta, lam, theta = estimates[:, -1]
    assert eta + (1 - eta) * math.exp(-lam) == pytest.approx((x == 0).mean(), abs=0.02)
    assert 0 <= eta <= 1 and lam > 0 and theta >= 0
    with pytest.raises(ValueError):
        e.zigp_em([3.0])   # one count: no variance for the starting values
    # Integers have R's mean of an integer vector, floats its mean of doubles.
    whole = np.array([0, 0, 1, 1, 2, 3, 5, 9])
    np.testing.assert_allclose(e.zigp_em(whole), e.zigp_em(whole + 0.0), rtol=1e-6)


def test_zigp_negloglik():
    count = np.array([5.0, 3.0, 2.0, 1.0])
    expected = -(count[:3] * np.log(dzigp(np.arange(3), 0.4, 2.0, 0.3))).sum()
    assert e.zigp_negloglik(count, 0.4, 2.0, 0.3, 2) == pytest.approx(expected, rel=1e-13)


def test_eb_fit_and_local_fdr():
    fit = e.eb_fit(COUNTS)
    K = int(COUNTS.max())
    assert fit["estimates"].shape == (3, K) and len(fit["pi0"]) == K and len(fit["criterion"]) == K
    assert 1 <= fit["cutoff"] <= K and (fit["pi0"] <= 1).all()
    assert fit["eta"] == fit["estimates"][0, fit["cutoff"] - 1]
    assert fit["pi0_hat"] == fit["pi0"][fit["cutoff"] - 1]
    assert 3 <= fit["cutoff"] <= 13   # between the two groups of counts
    fdr = e.eb_local_fdr(fit, [0.0, 20.0], [2.0, 70.0])
    f0 = dzigp([0, 20], fit["eta"], fit["lambda"], fit["theta"])
    np.testing.assert_allclose(fdr, fit["pi0_hat"] * (f0 + (1 - fit["eta"]) * 1e-10 + 1e-10) *
                               np.array([2.0, 70.0]), rtol=1e-12)
    assert e.eb_rejections(COUNTS) == 10.0   # the ten well-separated counts
    assert e.eb_rejections(COUNTS, level=0.0) == 0.0


def test_second_stage_threshold():
    fit = {"lambda": 2.0, "theta": 0.3, "cutoff": 5}
    bound = max(2.0 / (math.exp(0.3 - 1) - 0.3), math.log(100) / (0.3 - 1 - math.log(0.3)))
    assert e.second_stage_threshold(fit, 100) == math.ceil(max(6, bound))
    assert e.second_stage_threshold({"lambda": 0.1, "theta": 0.01, "cutoff": 9}, 10) == 10.0


def test_storey_and_simulation_counts():
    data = simulated_data(SETTING, 889192)
    fit = e.eb_fit(data["x"])
    rejected = e.storey_rejected(data, fit, 0.05)
    assert len(set(rejected)) == len(rejected) and (rejected < len(data["values"])).all()
    assert len(e.storey_rejected(data, fit, 0.0)) == 0
    assert (data["values"][rejected] > 5).mean() > 0.5
    counts = e.eb_simulation_counts(data, 5, 0.05)
    assert list(counts) == ["TS_R", "TS_V", "TS_TP", "OS_R", "OS_V", "OS_TP", "Storey_R",
                            "Storey_V", "Storey_TP"]
    for name in ("TS", "OS", "Storey"):
        assert counts[name + "_R"] == counts[name + "_V"] + counts[name + "_TP"]
        assert 0 <= counts[name + "_TP"] <= len(data["alt"])
    assert counts["Storey_V"] == data["null_frequency"][rejected].sum()
    assert counts["TS_TP"] >= 0.5 * len(data["alt"])


@pytest.mark.against_r
def test_same_fit_and_counts_as_r(run_r, reader):
    r = reader(run_r(r"""
setting <- function(N, lambda, theta, C) data.frame(N = N, eta = 0.4, lambda = lambda, theta = theta, pi0 = 0.8, f1 = "binomial", C = C)
cases <- list(list(setting(50, 2, 0.3, 5), 889192), list(setting(200, 4, 0.15, 5), 889200),
              list(setting(1000, 2, 0.3, 10), 889195), list(setting(500, 1.5, 0, 10), 3),
              list(setting(100, 1, 0.3, 5), 77))
for (case in cases) {
  data <- simulated_data(case[[1]], case[[2]]); fit <- eb_fit(data$x)
  put(fit$estimates); put(fit$pi0); put(fit$criterion); put(fit$cutoff)
  put(c(fit$eta, fit$lambda, fit$theta, fit$pi0_hat, second_stage_threshold(fit, length(data$x))))
  put(length(storey_rejected(data, fit, 0.05))); put(storey_rejected(data, fit, 0.05))
  put(unlist(eb_simulation_counts(data, case[[1]]$C, 0.05))); put(eb_rejections(data$x, 0.05))
  put(unlist(em_tail(max(data$x), 0.4, 2, 0.3, 3, 40))); put(zigp_negloglik(tabulate(data$x + 1), 0.4, 2, 0.3, 4))
}
for (domain in DOMAINS) { x <- read_domain(domain, PATHS$domains)$counts; fit <- eb_fit(x)
  put(fit$estimates); put(fit$criterion); put(fit$cutoff); put(eb_rejections(x, 0.05)) }
put(gp_pmf(0:60, 2.5, 0.35)); put(zigp_pmf_floored(0:60, 0.3, 2.5, 0.35)); put(zigp_cdf(0.3, 2.5, 0.35, 7))
"""))
    import os
    from localfdr.real_data import DOMAINS, read_domain
    from localfdr.setup import PATHS

    def setting(N, lam, theta, C):
        return {"N": N, "eta": 0.4, "lambda": lam, "theta": theta, "pi0": 0.8, "f1": "binomial", "C": C}

    cases = [(setting(50, 2, 0.3, 5), 889192), (setting(200, 4, 0.15, 5), 889200),
             (setting(1000, 2, 0.3, 10), 889195), (setting(500, 1.5, 0, 10), 3),
             (setting(100, 1, 0.3, 5), 77)]
    for case, seed in cases:
        data = simulated_data(case, seed)
        fit = e.eb_fit(data["x"])
        label = f"N = {case['N']}, seed {seed}: "
        r.same(label + "estimates", fit["estimates"])
        r.same(label + "pi0", fit["pi0"])
        r.same(label + "criterion", fit["criterion"])
        r.same(label + "cutoff", [fit["cutoff"]])
        r.same(label + "chosen", [fit["eta"], fit["lambda"], fit["theta"], fit["pi0_hat"],
                                  e.second_stage_threshold(fit, len(data["x"]))])
        storey = e.storey_rejected(data, fit, 0.05)
        r.same(label + "Storey count", [len(storey)])
        r.same(label + "Storey positions", storey + 1)
        r.same(label + "counts", list(e.eb_simulation_counts(data, case["C"], 0.05).values()))
        r.same(label + "rejections", [e.eb_rejections(data["x"], 0.05)])
        r.same(label + "em_tail", list(e.em_tail(data["x"].max(), 0.4, 2, 0.3, 3, 40).values()))
        count = np.bincount(data["x"].astype(int)).astype(float)
        r.same(label + "negloglik", [e.zigp_negloglik(count, 0.4, 2, 0.3, 4)])
    for domain in DOMAINS:
        x = read_domain(domain, PATHS["domains"])["counts"]
        fit = e.eb_fit(x)
        r.same(domain + " estimates", fit["estimates"])
        r.same(domain + " criterion", fit["criterion"])
        r.same(domain + " cutoff", [fit["cutoff"]])
        r.same(domain + " rejections", [e.eb_rejections(x, 0.05)])
    r.same("gp_pmf", e.gp_pmf(np.arange(61), 2.5, 0.35))
    r.same("zigp_pmf_floored", e.zigp_pmf_floored(np.arange(61), 0.3, 2.5, 0.35))
    r.same("zigp_cdf", [e.zigp_cdf(0.3, 2.5, 0.35, 7)])
    r.done()
    assert os.path.isdir(PATHS["domains"])
