# -----------------------------------------------------------------------------
# File:    python/tests/test_bayes_fdr.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Tests of bayes_fdr.py: every step of the sampler on small inputs with
#          a known answer, and seeded chains that give the same draws as R,
#          for every model and version of the sampler.
# Usage:   python -m pytest python/tests/test_bayes_fdr.py
# -----------------------------------------------------------------------------
import math

import numpy as np
import pytest

from inputs import COUNTS, state
from localfdr import rbase, rng
from localfdr import bayes_fdr as b
from localfdr.distributions import dzigp, shifted_density

V18, V07, V19 = (b.SAMPLER_VERSIONS[name] for name in ("2018-03-25", "2021-02-07", "2021-02-19"))
U_VALUES = np.unique(COUNTS)
FREQUENCY = np.array([(COUNTS == value).sum() for value in U_VALUES], dtype=float)


def test_sampler_settings_and_versions():
    settings = b.sampler_settings(1, 0.5)
    assert settings["init"]["lambda"] == 1 and settings["init"]["theta"] == 0.5
    assert len(settings["init"]["beta"]) == 50 and settings["init"]["beta"].sum() == pytest.approx(1)
    assert b.sampler_settings()["init"]["theta"] == 0.25
    assert settings["proposal"]["null_var"].shape == (3, 3)
    assert set(b.SAMPLER_VERSIONS) == {"2018-03-25", "2021-02-07", "2021-02-19"}
    assert V18["extra_atoms"](10, 40) == 60 and V19["extra_atoms"](10, 40) == 20
    assert V18["fdr_from"](1000) == 1 and V19["fdr_from"](1000) == 1001
    assert b._model_code("NP") == 3
    with pytest.raises(ValueError):
        b._model_code("XX")


def test_mixture_fdr():
    u = np.array([0.0, 1.0, 2.0, 7.0, 8.0])
    s = state(5, psi_np=np.array([0.1, 0.1, 0.1, 0.2, 0.0]), psi_sp=np.full(5, 0.3))
    f0 = 0.5 * dzigp(u, 0.5, 2.0, 0.25)
    out = b.mixture_fdr("NP", u, s, 2, V19)
    np.testing.assert_array_equal(out["f0"], f0)
    np.testing.assert_array_equal(out["f"], [f0[0], f0[1], f0[2], 0.2, 0.0])
    np.testing.assert_array_equal(out["fdr"], [1, 1, 1, min(1, f0[3] / 0.2), 1])
    out = b.mixture_fdr("SP", u, s, 2, V19)
    np.testing.assert_allclose(out["f"][3:], f0[3:] + 0.5 * 0.3, rtol=1e-15)
    out = b.mixture_fdr("P", u, s, 2, V19)
    np.testing.assert_allclose(out["f"], f0 + 0.5 * shifted_density("P", u, 1.0, 1.0, 0.5, 2))
    out = b.mixture_fdr("GP", u, s, None, V19)   # no cutoff: C = 1
    np.testing.assert_allclose(out["f"], f0 + 0.5 * shifted_density("GP", u, 1.0, 1.0, 0.5, 1))
    with pytest.raises(ValueError):
        b.mixture_fdr("GP", u, s, None, V18)
    short = b.mixture_fdr("NP", u, state(5, psi_np=np.full(4, 0.1)), 2, V19)
    assert np.isnan(short["f"][4]) and np.isnan(short["fdr"][4])


def test_draw_cutoff_and_tau():
    rng.set_seed(1)
    s = state(len(U_VALUES))
    draws = [b.draw_cutoff("NP", U_VALUES, FREQUENCY, s, {"kappa": 0.5, "zeta": 1e-10}, V19)
             for _ in range(50)]
    assert set(draws) <= set(U_VALUES)
    assert b.draw_cutoff("NP", U_VALUES, FREQUENCY, state(len(U_VALUES), tau=math.nan),
                         {"kappa": 0.5, "zeta": 1e-10}, V19) is None
    with pytest.raises(ValueError):
        b.draw_cutoff("NP", U_VALUES, FREQUENCY, state(3, psi_np=np.array([0.1, np.nan, 0.1])),
                      {"kappa": 0.5, "zeta": 1e-10}, V19)
    taus = np.array([b.draw_tau(7, {"kappa": 0.5, "zeta": 1e-10}) for _ in range(4000)])
    assert abs(taus.mean() - 7.5) < 0.2
    assert b.draw_tau(None, {"kappa": 0.5, "zeta": 1e-10}) > 0


def test_draw_groups_and_xi():
    rng.set_seed(2)
    x_sorted = np.sort(COUNTS)
    groups = b.draw_groups("NP", U_VALUES, x_sorted, FREQUENCY, state(len(U_VALUES)), 5, V19)
    assert groups["K"] == 20
    np.testing.assert_array_equal(groups["n0"] + groups["n1"], FREQUENCY)
    assert groups["N0"] == groups["n0"].sum() and groups["n00"] == groups["n0"][0] == 30
    assert groups["N1"] == groups["n1"][U_VALUES > 5].sum()
    assert groups["N2"] == (U_VALUES * groups["n0"]).sum()
    assert groups["N3"] == (U_VALUES * groups["n1"])[U_VALUES > 5].sum()
    assert (groups["n1"][U_VALUES <= 5] == 0).all()   # null with probability 1 up to the cutoff
    missing = b.draw_groups("NP", U_VALUES, x_sorted, FREQUENCY,
                            state(len(U_VALUES), pi0=math.nan), 5, V07)
    assert math.isnan(missing["N0"]) and np.isnan(missing["n1"]).all()
    replaced = b.draw_groups("NP", U_VALUES, x_sorted, FREQUENCY,
                             state(len(U_VALUES), pi0=math.nan), 5, V19)
    assert replaced["N0"] + replaced["n1"].sum() == len(COUNTS)
    groups = {"K": 20.0, "N3": 170.0, "N1": 10.0}
    draws = np.array([b.draw_xi(5, groups, {"kappa": 0.5, "zeta": 1e-10}) for _ in range(4000)])
    assert abs(draws.mean() - 111.5 / 10) < 0.1   # shape 170 + 0.5 - 6 * 10 + 1, rate 10


def test_concentration_above_cutoff_and_draw_atoms():
    np.testing.assert_array_equal(b.concentration_above_cutoff(2, 6, 1.5, 0.01),
                                  [0.01, 0.01, 0.01, 1.5, 1.5, 1.5])
    np.testing.assert_array_equal(b.concentration_above_cutoff(-1, 3, 1.5), [1.5, 1.5, 1.5])
    np.testing.assert_array_equal(b.concentration_above_cutoff(9, 3, 1.5), [1e-6, 1e-6, 1e-6])
    with pytest.raises(ValueError):
        b.concentration_above_cutoff(None, 3, 1.5)
    rng.set_seed(3)
    U = len(U_VALUES)
    s = state(U)
    psi = b.draw_atoms(U_VALUES, s["psi_np"], np.full(U, 1.0) + FREQUENCY, s, 1.0, 5, 30, V19)
    assert len(psi) == U
    np.testing.assert_array_equal(psi[:6], 0.5 * dzigp(U_VALUES[:6], 0.5, 2.0, 0.25))
    assert 0 < psi[6:].sum() <= 1 - 6 / 50 + 1e-12
    with pytest.raises(ValueError):
        b.draw_atoms(U_VALUES, s["psi_np"], np.full(U, 1.0), s, 1.0, None, 30, V18)
    assert len(b.draw_atoms(U_VALUES, s["psi_np"], np.full(U, 1.0), s, 1.0, None, 5, V19)) == U


def test_parameters_and_likelihoods():
    np.testing.assert_allclose(b.null_parameters([0.0, math.log(3.0), 0.0], V19), [0.5, 3.0, 0.5])
    assert b.null_parameters([0.0, 0.0, 1.0], V19)[2] == pytest.approx(1 / (1 + math.exp(-2)))
    np.testing.assert_allclose(b.gp_parameters([math.log(2.0), 0.0], V18), [2.0, 0.5])
    groups = {"n00": 3.0, "N0": 8.0, "N2": 9.0, "n0": np.array([3.0, 2.0, 2.0, 1.0]),
              "N1": 4.0, "N3": 11.0, "n1": np.array([0.0, 1.0, 1.0, 2.0])}
    eta, lam, theta = 0.3, 1.7, 0.2
    y = np.array([1.0, 2.0, 3.0])
    expected = 3 * math.log(eta + (1 - eta) * math.exp(-lam)) + 5 * math.log(1 - eta) + \
        5 * (math.log(lam) - lam) - 9 * theta + \
        (groups["n0"][1:] * (y - 1) * np.log(lam + theta * y)).sum() - \
        (groups["n0"][1:] * np.array([0.0, math.log(2), math.log(6)])).sum()
    assert b.loglik_null([eta, lam, theta], groups) == pytest.approx(expected, rel=1e-14)
    expected_gp = 4 * math.log(1.2) - 1.2 * 4 - 0.1 * 11 + \
        (groups["n1"][1:] * (y - 1) * np.log(1.2 + 0.1 * y)).sum() - \
        (groups["n1"][1:] * np.array([0.0, math.log(2), math.log(6)])).sum()
    assert b.loglik_gp([1.2, 0.1], groups) == pytest.approx(expected_gp, rel=1e-14)
    pre = [0.2, 0.1, -0.3]
    par = b.null_parameters(pre, V19)
    jacobian = math.exp(-0.2) / (1 + math.exp(-0.2)) ** 2 * math.exp(0.1) * \
        math.exp(-0.6) / (0.5 * (1 + math.exp(-0.6)) ** 2)
    expected = b.loglik_null(par, groups) + 0 - 0.5 * math.log(par[1]) + math.log(2) + \
        math.log(jacobian)
    assert b.logpost_null(pre, groups, V19) == pytest.approx(expected, rel=1e-13)
    assert b.logpost_null([0.2, 0.1, 0.3], groups, V19) == -math.inf   # theta above 0.5
    par = b.gp_parameters([0.2, -0.3], V18)
    expected = b.loglik_gp(par, groups) - 0.5 * math.log(par[0]) + 0.2 + \
        math.log(math.exp(-0.3) / (1 + math.exp(-0.3)) ** 2)
    assert b.logpost_gp([0.2, -0.3], groups, V18) == pytest.approx(expected, rel=1e-13)
    with pytest.raises(ValueError):
        b.loglik_null([eta, lam, theta], {**groups, "n0": np.array([3.0])})


def test_positive_definite_and_adapted_proposal():
    M = np.array([[2.0, 1.0], [1.0, 2.0]])
    np.testing.assert_allclose(b.positive_definite(M), M, rtol=1e-15)
    singular = np.array([[1.0, 1.0], [1.0, 1.0]])
    lifted = b.positive_definite(singular)
    assert np.linalg.eigvalsh(lifted).min() > 0 and np.abs(lifted - singular).max() < 1e-14
    with pytest.raises(ValueError):
        b.positive_definite(np.ones((2, 3)))
    rng.set_seed(4)
    history = rng.rnorm(600).reshape(200, 3)
    half = b.adapted_proposal(history, 200, {"prop": 0.5}, np.diag([0.1] * 3))
    np.testing.assert_allclose(half.T @ half, 100 / 101 * np.cov(history[99:200], rowvar=False),
                               rtol=1e-12)
    constant = np.ones((200, 3))
    fallback = b.adapted_proposal(constant, 200, {"prop": 0.5}, np.diag([0.1] * 3))
    assert np.isfinite(fallback).all()
    np.testing.assert_array_equal(b._proposal(np.zeros((2, 2)), np.array([1.0, 2.0])), [1.0, 2.0])
    assert b._counts({"n0": [1.0, 2.0]}, "n0").tolist() == [1.0, 2.0]
    assert b._vector(3).tolist() == [3.0]
    assert b._null_density([0.0], state(1), V19)[0] == 0.5 * dzigp(0, 0.5, 2.0, 0.25)[0]


@pytest.mark.parametrize("model", ["NP", "SP", "P", "GP"])
def test_bayes_local_fdr_separates_the_groups(model):
    rng.set_seed(10)
    fit = b.bayes_local_fdr(COUNTS, model, draws=400, burn=200)
    np.testing.assert_array_equal(fit["values"], U_VALUES)
    np.testing.assert_array_equal(fit["frequency"], FREQUENCY)
    assert (fit["fdr"][U_VALUES <= 3] > 0.9).all() and (fit["fdr"][U_VALUES >= 15] < 0.05).all()
    assert len(fit["eta"]) == 400 and len(fit["acceptance"]) == 600
    assert (fit["pi0"] > 0.6).all() and (fit["pi0"] < 1).all()
    rng.set_seed(10)
    again = b.bayes_local_fdr(COUNTS, model, draws=400, burn=200)
    np.testing.assert_array_equal(again["fdr"], fit["fdr"])


def test_bayes_local_fdr_options():
    rng.set_seed(11)
    fit = b.bayes_local_fdr(COUNTS, "NP", draws=60, burn=20, zigp=False, cutoff=False, thin=2)
    assert (fit["theta"] == 0).all() and (fit["C"] == -1).all() and len(fit["acceptance"]) == 140
    with pytest.raises(IndexError):
        b.bayes_local_fdr(COUNTS, "NP", draws=10, burn=20)   # no stored draw after the burn-in
    with pytest.raises(ValueError):
        b.bayes_local_fdr(COUNTS, "XX")


R_CHAINS = r"""
setting <- function(N, C = 5) data.frame(N = N, eta = 0.4, lambda = 2, theta = 0.3, pi0 = 0.8, f1 = "binomial", C = C)
run <- function(x, model, version, seed, ...) {
  set.seed(seed)
  fit <- tryCatch(suppressWarnings(bayes_local_fdr(x, model, draws = 260, burn = 140, version = version, ...)),
                  error = function(e) NULL)
  if (is.null(fit)) return(put(-999))
  put(1)
  for (name in c("values", "frequency", "fdr", "eta", "lambda", "theta", "pi0", "tau", "C", "xi", "delta", "nu", "acceptance")) put(fit[[name]])
  put(runif(1))
}
x1 <- simulated_data(setting(100), 7)$x
x2 <- simulated_data(setting(1000, 10), 1441)$x
x3 <- read_domain("cd00083", "data/domains")$counts
x4 <- c(rep(0:3, c(40, 20, 10, 5)), 10:65)
for (version in names(SAMPLER_VERSIONS)) for (model in c("NP", "SP", "P", "GP")) {
  run(x1, model, version, 11, gamma = 1)
  run(x3, model, version, 13, gamma = 1.5, settings = sampler_settings(1, 0.5))
  run(x4, model, version, 14, gamma = 0.5)
}
run(x2, "NP", "2021-02-19", 12); run(x2, "GP", "2021-02-07", 12)
run(x1, "NP", "2021-02-19", 5, cutoff = FALSE); run(x1, "GP", "2021-02-19", 5, zigp = FALSE)
run(x1, "SP", "2021-02-07", 5, thin = 2, gamma = 0.02)
"""
NAMES = ("values", "frequency", "fdr", "eta", "lambda", "theta", "pi0", "tau", "C", "xi", "delta",
         "nu", "acceptance")


@pytest.mark.against_r
def test_seeded_chains_give_the_draws_of_r(run_r, reader):
    """Every stored draw of every parameter, the acceptance rate of every
    iteration, the posterior mean fdr and the next number of the stream:
    three versions x four models x three data sets (simulated, a protein
    domain, 60 distinct counts), and the options. Where R stops with an
    error, so does Python."""
    import os
    from localfdr.real_data import read_domain
    from localfdr.setup import REPO_ROOT
    from localfdr.simulate import simulated_data

    r = reader(run_r(R_CHAINS))

    def setting(N, C=5):
        return {"N": N, "eta": 0.4, "lambda": 2, "theta": 0.3, "pi0": 0.8, "f1": "binomial", "C": C}

    def run(label, x, model, version, seed, **options):
        rng.set_seed(seed)
        try:
            fit = b.bayes_local_fdr(x, model, draws=260, burn=140, version=version, **options)
        except (ValueError, IndexError):
            fit = None
        flag = r.take(1)[0]
        assert (flag == -999) == (fit is None), f"{label}: an error in one language only"
        if fit is None:
            return "error"
        for name in NAMES:
            r.same(f"{label}: {name}", fit[name])
        r.same(f"{label}: stream after the chain", rng.runif(1))
        return "missing" if np.isnan(fit["fdr"]).any() else "ok"

    x1 = simulated_data(setting(100), 7)["x"]
    x2 = simulated_data(setting(1000, 10), 1441)["x"]
    x3 = read_domain("cd00083", os.path.join(REPO_ROOT, "data", "domains"))["counts"]
    x4 = np.concatenate((np.repeat(np.arange(4.0), [40, 20, 10, 5]), np.arange(10.0, 66.0)))
    outcomes = []
    for version in b.SAMPLER_VERSIONS:
        for model in ("NP", "SP", "P", "GP"):
            outcomes.append(run(f"{version} {model} simulated", x1, model, version, 11, gamma=1))
            outcomes.append(run(f"{version} {model} domain", x3, model, version, 13, gamma=1.5,
                                settings=b.sampler_settings(1, 0.5)))
            outcomes.append(run(f"{version} {model} 60 counts", x4, model, version, 14, gamma=0.5))
    outcomes.append(run("N = 1000 NP", x2, "NP", "2021-02-19", 12))
    outcomes.append(run("N = 1000 GP", x2, "GP", "2021-02-07", 12))
    outcomes.append(run("no cutoff", x1, "NP", "2021-02-19", 5, cutoff=False))
    outcomes.append(run("zero-inflated Poisson", x1, "GP", "2021-02-19", 5, zigp=False))
    outcomes.append(run("thinned", x1, "SP", "2021-02-07", 5, thin=2, gamma=0.02))
    r.done()
    assert outcomes.count("ok") >= 30
