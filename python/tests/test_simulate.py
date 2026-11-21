# -----------------------------------------------------------------------------
# File:    python/tests/test_simulate.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Tests of simulate.py: the simulated data, the settings of the runs,
#          one replication of each kind, the driver and the summaries on
#          hand-made replications, the layout of the tables, and the same
#          data, replications and summaries as R.
# Usage:   python -m pytest python/tests/test_simulate.py
# -----------------------------------------------------------------------------
import os

import numpy as np
import pandas as pd
import pytest

from inputs import SETTING, SHORT
from localfdr import rbase, rng
from localfdr import simulate as sim
from localfdr.setup import PAPER, PATHS


def test_rmixture():
    data = sim.rmixture(400, 0.4, 2.0, 0.3, 0.8, "binomial", 5, seed=3)
    assert len(data["draw"]) == 400 and set(data["true_grp"]) == {0.0, 1.0}
    assert (data["draw"][data["true_grp"] == 1] > 5).all()
    assert abs(data["true_grp"].mean() - 0.2) < 0.06
    again = sim.rmixture(400, 0.4, 2.0, 0.3, 0.8, "binomial", 5, seed=3)
    np.testing.assert_array_equal(again["draw"], data["draw"])
    for f1 in ("geometric", "binomial mixture", "poisson mixture"):
        other = sim.rmixture(300, 0.4, 2.0, 0.3, 0.5, f1, 5, seed=4)
        assert (other["draw"][other["true_grp"] == 1] > 5).all()
    assert (sim.rmixture(50, 0.4, 2.0, 0.3, 1.0, "unknown", 5, seed=1)["true_grp"] == 0).all()
    with pytest.raises(ValueError):
        sim.rmixture(50, 0.4, 2.0, 0.3, 0.5, "unknown", 5, seed=1)
    with pytest.raises(ValueError):
        sim.rmixture(50, 0.4, 2.0, 0.3, 1.5, "binomial", 5, seed=1)


def test_simulated_data():
    data = sim.simulated_data(SETTING, 889192)
    assert len(data["x"]) == 100 and len(data["null"]) + len(data["alt"]) == 100
    np.testing.assert_array_equal(data["x"], np.concatenate((data["null"], data["alt"])))
    np.testing.assert_array_equal(data["values"], np.unique(data["x"]))
    assert data["null_frequency"].sum() == len(data["null"])
    assert data["alt_frequency"].sum() == len(data["alt"])
    assert (data["alt_frequency"][data["values"] <= 5] == 0).all()
    with pytest.raises(ValueError):
        sim.simulated_data({**SETTING, "N": 3, "pi0": 1.0}, 1)


def test_settings_of_the_runs():
    every = sim.simulation_settings()
    assert list(every.columns) == ["run", "setting", "eta", "lambda", "theta", "pi0", "f1", "C", "N"]
    assert len(every) == 93 and not every.duplicated(["run", "setting", "N"]).any()
    assert (every["eta"] == 0.4).all() and set(every["pi0"]) == {0.8, 0.35}
    assert sim._grid("r", [1, 2], 3, [0.1, 0.2], 5, 50).shape == (2, 9)
    row = every[(every["run"] == "run0219") & (every["setting"] == 17)].iloc[0]
    assert (row["lambda"], row["theta"], row["N"], row["C"]) == (2, 0.3, 1000, 5)
    assert len(sim.item_settings("table1")) == 10 and len(sim.item_settings("table2")) == 5
    assert len(sim.item_settings("figS1")) == 22 and len(sim.item_settings("figS3")) == 60
    assert set(sim.item_settings("table2")["C"]) == {10}
    with pytest.raises(ValueError):
        sim.item_settings("table9")
    assert set(sim.SIMULATION_RUNS) == set(every["run"])
    assert all(run["version"] in sim.SAMPLER_VERSIONS for run in sim.SIMULATION_RUNS.values())


def test_replications_and_their_folder():
    out = sim.bayes_replication(SETTING, "NP", 1, 3, SHORT)
    assert list(out) == ["d", "R", "V", "TP", "n_alt"] and out["d"] == 3
    assert out["R"] == out["V"] + out["TP"] and 0 <= out["TP"] <= out["n_alt"]
    assert out == sim.bayes_replication(SETTING, "NP", 1, 3, SHORT)   # the seed is that of the data
    eb = sim.eb_replication(SETTING, 3, SHORT)
    assert eb["n_alt"] == out["n_alt"] and eb["TS_R"] == eb["TS_V"] + eb["TS_TP"]
    assert sim.replication_dir(SETTING, "EB", np.nan, PAPER, "out") == os.path.join(
        "out", "replications", "eta0.4_lambda2_theta0.3_pi00.8_binomial_C5_N100", "EB_seed889191")
    assert sim.replication_dir(SETTING, "NP", 1.5, PAPER, "out").endswith(
        "NP_gamma1.5_seed889191_sampler2021-02-19_init0.25_draws5000_burn1000")
    assert sim._replication_file({**SETTING, "method": "EB", "gamma": np.nan}, 7, PAPER,
                                 "out").endswith(os.path.join("EB_seed889191", "rep0007.csv"))


def test_driver_and_summaries(tmp_path):
    cells = pd.DataFrame([{**SETTING, "method": "NP", "gamma": 1.0},
                          {**SETTING, "method": "EB", "gamma": np.nan}])
    sim.run_replications(cells, [1, 2], SHORT, str(tmp_path), cores=2)
    folder = sim.replication_dir(SETTING, "NP", 1.0, SHORT, str(tmp_path))
    assert sorted(os.listdir(folder)) == ["rep0001.csv", "rep0002.csv"]
    first = open(os.path.join(folder, "rep0001.csv")).read()
    assert first.splitlines()[0] == '"d","R","V","TP","n_alt"'
    sim.run_replications(cells, [1, 2, 3], SHORT, str(tmp_path), cores=1)   # runs only the third
    assert open(os.path.join(folder, "rep0001.csv")).read() == first
    summary = sim.summarize_replications(cells, 3, SHORT, str(tmp_path))
    assert list(summary["method"]) == ["NP", "TS", "OS", "Storey"]
    reps = pd.concat([rbase.read_csv(os.path.join(folder, f"rep000{d}.csv")) for d in (1, 2, 3)])
    assert summary["R"][0] == round(reps["R"].mean() + 1e-12, 2)
    assert summary["TPR_sd"][0] == pytest.approx((reps["TP"] / reps["n_alt"]).std(), abs=5.1e-5)
    with pytest.raises(FileNotFoundError):
        sim.summarize_replications(cells, 4, SHORT, str(tmp_path))
    failing = pd.DataFrame([{**SETTING, "N": 3, "pi0": 1.0, "method": "EB", "gamma": np.nan}])
    with pytest.raises(RuntimeError):
        sim.run_replications(failing, [1], SHORT, str(tmp_path), cores=1)


def test_summaries_of_hand_made_replications(tmp_path):
    cells = pd.DataFrame([{**SETTING, "method": "GP", "gamma": 0.5},
                          {**SETTING, "method": "EB", "gamma": np.nan}])
    rows = [(4, 1, 3, 10), (0, 0, 0, 8), (6, 0, 6, 12)]
    for d, (R, V, TP, n_alt) in enumerate(rows, start=1):
        for cell in cells.to_dict("records"):
            file = sim._replication_file(cell, d, PAPER, str(tmp_path))
            os.makedirs(os.path.dirname(file), exist_ok=True)
            row = {"d": d, "R": R, "V": V, "TP": TP, "n_alt": n_alt} if cell["method"] == "GP" else \
                {"d": d, "n_alt": n_alt, "TS_R": R, "TS_V": V, "TS_TP": TP, "OS_R": R, "OS_V": V,
                 "OS_TP": TP, "Storey_R": R, "Storey_V": V, "Storey_TP": TP}
            rbase.write_csv(pd.DataFrame([row]), file)
    summary = sim.summarize_replications(cells, 3, PAPER, str(tmp_path)).set_index("method")
    assert summary.loc["GP", "R"] == 3.33 and summary.loc["GP", "R_sd"] == 3.06
    assert summary.loc["GP", "FDR"] == 0.0833        # (1/4 + 0 + 0) / 3
    assert summary.loc["Storey", "FDR"] == 0.125     # the replication without rejections left out
    assert summary.loc["GP", "TPR"] == 0.2667 and summary.loc["TS", "TPR_sd"] == 0.2517
    assert np.isnan(summary.loc["OS", "gamma"]) and summary.loc["GP", "gamma"] == 0.5


def test_simulation_summaries(tmp_path, capsys):
    paths = {**PATHS, "results": str(tmp_path)}
    settings = {**SHORT, "sim_draws": 40, "sim_burn": 20}
    opts = {"runs": "1", "replications": 2, "cores": 2}
    assert sim.simulation_summaries("table2", opts, settings, paths) is None
    assert "not all stored yet" in capsys.readouterr().out
    summary = sim.simulation_summaries("table2", {**opts, "runs": "2"}, settings, paths)
    assert len(summary) == 50 and set(summary["method"]) == {"NP", "SP"}
    table = sim.table2_layout(summary)
    assert table.shape == (25, 15) and not table.isna().any().any()


def test_archive_and_layouts():
    archive = sim.archived_summaries("table1", PATHS)
    assert set(archive["study"]) == {"methods"} and set(archive["run"]) == {"run0219", "run0224"}
    assert len(archive) == 60
    with_c = sim.with_parameters(archive)
    assert set(with_c["C"]) == {5, 10} and (with_c["pi0"] == 0.8).all()
    block = sim.summary_block(archive.iloc[:3], "X")
    assert list(block.columns) == ["X_R", "X_R_sd", "X_FDR", "X_FDR_sd", "X_TPR", "X_TPR_sd"]
    table = sim.table1_layout(archive)
    assert table.shape == (30, 14) and list(table["N"][:7]) == [50] * 6 + [100]
    assert list(table["method"][:6]) == ["Nonparametric", "Parametric", "Semiparametric",
                                         "Two-Stage", "One-Stage", "Storey"]
    assert (table["C10_R"][0], table["C5_TPR_sd"][0]) == (9.33, 0.2104)   # Table 1, first row
    prior = sim.archived_summaries("table2", PATHS)
    table = sim.table2_layout(prior)
    assert table.shape == (25, 15) and list(table["prior"][:5]) == list(sim.PRIOR_GAMMAS)
    assert (table["NP_R"][0], table["SP_TPR_sd"][0]) == (10.7, 0.1906)    # Table 2, first row
    incomplete = sim.table1_layout(archive[archive["method"] != "GP"])
    assert incomplete["C10_R"].isna().sum() == 5
    assert sim._match([("a",), ("z",)], [("b",), ("a",), ("a",)]) == [1, -1]
    assert sim._take(archive, [0, -1])["R"].isna().tolist() == [False, True]
    assert np.isinf(sim._ratio(np.array([1.0]), np.array([0.0]))[0])


@pytest.mark.against_r
def test_same_data_replications_and_summaries_as_r(run_r, reader, tmp_path):
    r_results, py_results = tmp_path / "r", tmp_path / "py"
    r = reader(run_r(r"""
for (f1 in c("binomial", "geometric", "binomial mixture", "poisson mixture")) for (seed in c(1439, 889192, 5)) {
  data <- rmixture(600, 0.4, 2.5, 0.3, 0.7, f1, 5, seed); put(data$draw); put(data$true_grp); put(runif(1)) }
setting <- simulation_settings()[c(47, 23), ]
for (i in 1:2) { data <- simulated_data(setting[i, ], 1441); put(length(data$values))
  put(data$values); put(data$null_frequency); put(data$alt_frequency) }
settings <- PAPER; settings$sim_draws <- 130; settings$sim_burn <- 110
cells <- rbind(cbind(setting, method = "NP", gamma = 0.02), cbind(setting, method = "SP", gamma = 1),
               cbind(setting, method = "GP", gamma = 0.5), cbind(setting, method = "EB", gamma = NA))
run_replications(cells, 1:3, settings, "%s", cores = 2)
s <- summarize_replications(cells, 3, settings, "%s")
write.csv(s, "%s", row.names = FALSE)
""" % (r_results, r_results, tmp_path / "r_summary.csv")))
    for f1 in ("binomial", "geometric", "binomial mixture", "poisson mixture"):
        for seed in (1439, 889192, 5):
            data = sim.rmixture(600, 0.4, 2.5, 0.3, 0.7, f1, 5, seed)
            r.same(f"{f1} {seed}: counts", data["draw"])
            r.same(f"{f1} {seed}: groups", data["true_grp"])
            r.same(f"{f1} {seed}: stream", rng.runif(1))
    setting = sim.simulation_settings().iloc[[46, 22]]
    for row in setting.to_dict("records"):
        data = sim.simulated_data(row, 1441)
        r.same("number of values", [len(data["values"])])
        r.same("values", data["values"])
        r.same("null frequency", data["null_frequency"])
        r.same("non-null frequency", data["alt_frequency"])
    r.done()
    settings = {**PAPER, "sim_draws": 130, "sim_burn": 110}
    cells = pd.concat([setting.assign(method=method, gamma=gamma) for method, gamma in
                       (("NP", 0.02), ("SP", 1), ("GP", 0.5), ("EB", np.nan))], ignore_index=True)
    sim.run_replications(cells, [1, 2, 3], settings, str(py_results), cores=2)
    files = sorted(os.path.relpath(os.path.join(folder, name), r_results)
                   for folder, _, names in os.walk(r_results) for name in names)
    assert len(files) == 24
    for file in files:   # every replication file, byte for byte
        assert open(r_results / file).read() == open(py_results / file).read(), file
    rbase.write_csv(sim.summarize_replications(cells, 3, settings, str(py_results)),
                    tmp_path / "py_summary.csv")
    assert open(tmp_path / "py_summary.csv").read() == open(tmp_path / "r_summary.csv").read()
