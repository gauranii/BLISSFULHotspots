# -----------------------------------------------------------------------------
# File:    python/localfdr/simulate.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Tables 1 and 2, Figures S1 to S3.
# Purpose: The simulation study: data from the two-group mixture, the
#          settings of the four simulation runs behind the paper, one
#          replication of a Bayesian or an empirical Bayes method, the driver
#          that runs and stores replications, and the summaries (R, FDR, TPR)
#          from the archive or from stored replications.
# Inputs:  None.
# Outputs: None.
# Usage:   from localfdr.simulate import simulated_data, simulation_settings
# Source:  R/simulate.R. A seed gives the same data set as in R (rng.py), and
#          the summaries are rounded and written as R rounds and writes them.
# -----------------------------------------------------------------------------
import os

import numpy as np
import pandas as pd
from numba import njit

from . import rbase, rmath, rng
from .bayes_fdr import SAMPLER_VERSIONS, bayes_local_fdr, sampler_settings
from .distributions import _rzigp, check_zigp
from .eb_fdr import eb_simulation_counts
from .rng import _rbinom, _rgeom, _rpois, _unif
from .setup import log_message, parse_index_set

# ----- Data -----

F1_KINDS = ("binomial", "geometric", "binomial mixture", "poisson mixture")


@njit(cache=True, error_model="numpy")
def _rmixture(state, n, eta, lam, theta, pi0, kind, C, geometric_prob, draw, true_grp):
    """Fill draw and true_grp; returns 0, 1 when a ZIGP probability is
    missing, 2 when a non-null count is asked from an unknown f1."""
    first_component = 0.6  # weight of the first component of a mixture f1
    j = 0
    while j < n:
        if _unif(state) >= pi0:
            keep = True
            if kind == 0:
                zigp = _rzigp(state, eta, lam, theta)
                count = zigp + _rbinom(state, 250.0, 0.20) + C - 39
                keep = count > C
            elif kind == 1:
                zigp = _rzigp(state, eta, lam, theta)
                count = zigp + _rgeom(state, geometric_prob) + C + 1
            elif kind == 2:
                if _unif(state) >= first_component:
                    zigp = _rzigp(state, eta, lam, theta)
                    count = zigp + _rbinom(state, 250.0, 0.2) + C + 1 - 40
                    keep = count > C
                else:
                    zigp = _rzigp(state, eta, lam, theta)
                    count = zigp + _rbinom(state, 50.0, 0.2) + C + 1 + 0
            elif kind == 3:
                mean = 5.0 if _unif(state) >= first_component else 10.0
                zigp = _rzigp(state, eta, lam, theta)
                count = zigp + _rpois(state, mean) + C + 1
            else:
                return 2
            if zigp < 0:
                return 1
            if keep:
                draw[j] = count
                true_grp[j] = 1.0
                j += 1
        else:
            zigp = _rzigp(state, eta, lam, theta)
            if zigp < 0:
                return 1
            draw[j] = zigp
            true_grp[j] = 0.0
            j += 1
    return 0


def rmixture(n, eta, lam, theta, pi0, f1="binomial", C=None, seed=None, geometric_prob=0.25):
    """Simulate counts from the two-group mixture pi0 f0 + (1 - pi0) f1.

    A null count is a ZIGP draw. A non-null count is a ZIGP draw plus a
    shifted draw from f1; with a binomial f1 a non-null count at or below
    the cutoff is drawn again.

    Parameters: n (number of counts), eta, lam, theta (parameters of the
    ZIGP null f0), pi0 (proportion of null counts), f1 (alternative:
    "binomial", Binomial(250, 0.2) + C - 39; "geometric", Geometric + C + 1;
    "binomial mixture" or "poisson mixture"), C (cutoff: non-null counts lie
    above it), seed (seed, set at the start), geometric_prob (success
    probability of the geometric f1).
    Returns: dict with `draw` (the counts) and `true_grp` (1 for non-null).
    """
    rng.set_seed(seed)
    check_zigp(eta, lam, theta)
    if not (pi0 >= 0 and pi0 <= 1):
        raise ValueError("pi0 >= 0 && pi0 <= 1 is not TRUE")
    n = int(n)
    draw = np.zeros(n)
    true_grp = np.zeros(n)
    kind = F1_KINDS.index(f1) if f1 in F1_KINDS else -1
    status = _rmixture(rng._state, n, float(eta), float(lam), float(theta), float(pi0), kind,
                       float(C), float(geometric_prob), draw, true_grp)
    if status == 1:
        raise ValueError("missing value where TRUE/FALSE needed")
    if status == 2:
        raise ValueError(f"unknown f1: {f1}")
    return {"draw": draw, "true_grp": true_grp}


def simulated_data(setting, seed):
    """One simulated data set of a setting, split into null and non-null
    counts.

    Parameters: setting (one row of `simulation_settings()`), seed (seed of
    the data set).
    Returns: dict: `null` and `alt` (the counts of each group, in the order
    drawn), `x` (null counts followed by non-null counts) and `values`,
    `null_frequency`, `alt_frequency` (distinct counts and the number of
    null and non-null observations at each).
    """
    data = rmixture(setting["N"], setting["eta"], setting["lambda"], setting["theta"],
                    setting["pi0"], setting["f1"], setting["C"], seed)
    null = data["draw"][data["true_grp"] == 0]
    alt = data["draw"][data["true_grp"] == 1]
    if len(null) < 2 or len(alt) < 2:
        raise ValueError("a group has fewer than two counts")
    x = np.concatenate((null, alt))
    values = np.unique(x)

    def frequency(group):
        return np.bincount(np.searchsorted(values, group), minlength=len(values)).astype(np.float64)

    return {"null": null, "alt": alt, "x": x, "values": values,
            "null_frequency": frequency(null), "alt_frequency": frequency(alt)}


# ----- Settings of the simulation runs -----

# The four simulation runs behind the paper. `version` is the version of the
# sampler (SAMPLER_VERSIONS), `bayes_seed` and `eb_seed` the seeds of the data
# sets of the Bayesian and of the empirical Bayes methods (replication d uses
# seed + d), `lambda_init` and `theta_init` starting values of the sampler.
SIMULATION_RUNS = {
    "run0207": {"version": "2021-02-07", "bayes_seed": 889191, "eb_seed": 889191,
                "lambda_init": 2, "theta_init": 0.20},
    "run0219": {"version": "2021-02-19", "bayes_seed": 889191, "eb_seed": 889191,
                "lambda_init": 2, "theta_init": 0.25},
    "run0224": {"version": "2021-02-19", "bayes_seed": 1440, "eb_seed": 889191,
                "lambda_init": 2, "theta_init": 0.25},
    "run0227": {"version": "2021-02-19", "bayes_seed": 889191, "eb_seed": 889191,
                "lambda_init": 2, "theta_init": 0.25},
}

# Dirichlet concentrations of the prior sensitivity analysis (Table 2).
PRIOR_GAMMAS = {"Objective": 0.02, "Jeffreys'": 0.5, "Uniform": 1, "Near Uniform": 1.5,
                "Informative": 3}

# The six methods, as named in the archive, and the Dirichlet concentration
# (NP, SP) or unused prior argument (GP) of the Bayesian ones in Table 1.
METHODS = ("NP", "GP", "SP", "TS", "OS", "Storey")
METHOD_LABELS = {"NP": "Nonparametric", "GP": "Parametric", "SP": "Semiparametric",
                 "TS": "Two-Stage", "OS": "One-Stage", "Storey": "Storey"}
DEFAULT_GAMMA = {"NP": 1, "GP": 0.5, "SP": 1}
SUMMARY_COLUMNS = ("R", "R_sd", "FDR", "FDR_sd", "TPR", "TPR_sd")
_N_TABLE = (50, 100, 200, 500, 1000)


def _grid(run, setting, lam, theta, C, N, pi0=0.8):
    """Settings of one run that differ in the vectors given (recycled)."""
    length = max(np.size(setting), np.size(lam), np.size(theta), np.size(N))

    def column(value):
        return list(np.broadcast_to(np.asarray(value), (length,)))

    return pd.DataFrame({"run": run, "setting": column(setting), "eta": 0.4, "lambda": column(lam),
                         "theta": column(theta), "pi0": pi0, "f1": "binomial", "C": C,
                         "N": column(N)})


def simulation_settings():
    """Settings of the simulation runs: one row per run, setting number and N.

    Every setting has eta = 0.4, pi0 = 0.8 or 0.35, a shifted binomial f1
    and 1000 replications. The setting numbers are those of the settings
    file of each run.

    Returns: DataFrame with run, setting, eta, lambda, theta, pi0, f1, C, N.
    """
    # run0207: one setting, many N. Settings 16 and 18: Figure S1; 10 and 14: Figure S3.
    n_all = (50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
    parts = [_grid("run0207", 10, 4.0, 0.30, 5, _N_TABLE),
             _grid("run0207", 14, 1.0, 0.20, 5, _N_TABLE),
             _grid("run0207", 16, 1.5, 0.00, 10, n_all, pi0=0.35),
             _grid("run0207", 18, 1.5, 0.00, 10, n_all)]
    # run0219: settings 1 to 30, C = 5; six (lambda, theta) pairs for each N.
    lam, theta = (3, 4, 3, 4, 2, 3), (0.15, 0.15, 0.20, 0.20, 0.30, 0.30)
    for block, N in enumerate((100, 500, 1000, 50, 200)):
        parts.append(_grid("run0219", 6 * block + np.arange(1, 7), lam, theta, 5, N))
    # run0224: settings 1 to 11, C = 10, lambda = 2, theta = 0.3.
    parts.append(_grid("run0224", np.arange(1, 12), 2, 0.30, 10, n_all))
    # run0227: settings 1 to 20, C = 5; four (lambda, theta) pairs for each N.
    lam, theta = (1, 2, 2, 1), (0.15, 0.15, 0.20, 0.30)
    for block, N in enumerate(_N_TABLE):
        parts.append(_grid("run0227", 4 * block + np.arange(1, 5), lam, theta, 5, N))
    settings = pd.concat(parts, ignore_index=True)
    return settings.astype({"setting": np.int64, "lambda": np.float64, "theta": np.float64,
                            "C": np.int64, "N": np.int64})


def item_settings(item):
    """Rows of `simulation_settings()` behind one table or figure.

    Parameters: item ("table1", "table2", "figS1" or "figS3").
    Returns: DataFrame of settings.
    """
    every = simulation_settings()
    run, in_table = every["run"], every["N"].isin(_N_TABLE)
    if item == "table1":
        keep = ((run == "run0224") & in_table) | \
            ((run == "run0219") & (every["lambda"] == 2) & (every["theta"] == 0.3))
    elif item == "table2":
        keep = (run == "run0224") & in_table
    elif item == "figS1":
        keep = (run == "run0207") & every["setting"].isin((16, 18))
    elif item == "figS3":
        keep = ((run == "run0207") & every["setting"].isin((10, 14))) | \
            run.isin(("run0219", "run0227"))
    else:
        raise ValueError(f"unknown item: {item}")
    return every[keep]


# ----- One replication -----

def bayes_replication(setting, model, gamma, d, settings):
    """Rejections of a Bayesian method in one simulated data set.

    Parameters: setting (one row of `simulation_settings()`), model ("NP",
    "SP", "P" or "GP"), gamma (Dirichlet concentration), d (replication
    number; the data set has seed `bayes_seed + d`), settings (settings of
    the paper: `alpha`, `sim_draws`, `sim_burn`).
    Returns: dict (one row): d, R (rejections), V (false rejections), TP
    (true rejections), n_alt (non-null counts).
    """
    run = SIMULATION_RUNS[setting["run"]]
    data = simulated_data(setting, run["bayes_seed"] + d)
    # No second seed: as in the original runs the sampler continues the stream
    # that simulated_data() has just seeded, so seed + d fixes the data and the chain.
    fit = bayes_local_fdr(data["x"], model, gamma,
                          sampler_settings(run["lambda_init"], run["theta_init"]),
                          draws=settings["sim_draws"], burn=settings["sim_burn"],
                          version=run["version"])
    rejected = fit["fdr"] <= settings["alpha"]
    V = rbase.r_sum(data["null_frequency"][rejected])
    TP = rbase.r_sum(data["alt_frequency"][rejected])
    return {"d": int(d), "R": TP + V, "V": V, "TP": TP, "n_alt": len(data["alt"])}


def eb_replication(setting, d, settings):
    """Rejections of the three empirical Bayes methods in one simulated data
    set.

    Parameters: setting (one row of `simulation_settings()`), d (replication
    number; the data set has seed `eb_seed + d`), settings (settings of the
    paper: `alpha`).
    Returns: dict (one row): d, n_alt, and for each of TS (two-stage), OS
    (one-stage) and Storey the rejections `_R`, false rejections `_V` and
    true rejections `_TP`.
    """
    data = simulated_data(setting, SIMULATION_RUNS[setting["run"]]["eb_seed"] + d)
    return {"d": int(d), "n_alt": len(data["alt"]),
            **eb_simulation_counts(data, setting["C"], settings["alpha"])}


# ----- Driver -----

def replication_dir(setting, method, gamma, settings, results):
    """Folder of the stored replications of a setting and method.

    Settings of different runs that share data sets and sampler version
    share their replications.

    Parameters: setting (one row of `simulation_settings()`), method ("EB",
    or a Bayesian model), gamma (Dirichlet concentration, Bayesian methods),
    settings (settings of the paper), results (results folder).
    Returns: path of the folder.
    """
    run = SIMULATION_RUNS[setting["run"]]
    data = "eta%g_lambda%g_theta%g_pi0%g_%s_C%g_N%g" % (
        setting["eta"], setting["lambda"], setting["theta"], setting["pi0"], setting["f1"],
        setting["C"], setting["N"])
    if method == "EB":
        how = "EB_seed%d" % run["eb_seed"]
    else:
        how = "%s_gamma%g_seed%d_sampler%s_init%g_draws%d_burn%d" % (
            method, gamma, run["bayes_seed"], run["version"], run["theta_init"],
            settings["sim_draws"], settings["sim_burn"])
    return os.path.join(results, "replications", data, how)


def _replication_file(cell, d, settings, results):
    return os.path.join(replication_dir(cell, cell["method"], cell["gamma"], settings, results),
                        "rep%04d.csv" % d)


def run_replications(cells, runs, settings, results, cores=1):
    """Run the replications of a set of settings that are not yet stored.

    One csv file is written per setting, method and replication, so a study
    can be run in pieces and continued.

    Parameters: cells (DataFrame of settings with the columns `method`,
    "EB" or a Bayesian model, and `gamma`), runs (replication numbers to
    run), settings (settings of the paper), results (results folder), cores
    (number of parallel processes).
    Returns: None.
    """
    cells = cells.to_dict("records")
    jobs = [(cell, d, _replication_file(cell, d, settings, results)) for d in runs for cell in cells]
    jobs = [job for job in jobs if not os.path.exists(job[2])]
    log_message("%d replications to run, %d cores", len(jobs), cores)

    def run(job):
        cell, d, file = job
        if cell["method"] == "EB":
            result = eb_replication(cell, d, settings)
        else:
            result = bayes_replication(cell, cell["method"], cell["gamma"], d, settings)
        os.makedirs(os.path.dirname(file), exist_ok=True)
        rbase.write_csv(pd.DataFrame([result]), file)
        return True

    status = rbase.parallel_map(run, jobs, cores)
    failed = [result for result in status if result is not True]
    if failed:
        raise RuntimeError(f"{len(failed)} replications failed; first: {failed[0]!r}")


def _ratio(numerator, denominator):
    with np.errstate(divide="ignore", invalid="ignore"):
        return numerator / denominator


def summarize_replications(cells, replications, settings, results):
    """Mean and standard deviation of R, FDR and TPR over stored replications.

    Rounded as in the paper: R to two decimals, FDR and TPR to four. The FDR
    of a replication without rejections is 0 (Storey: left out of the mean).

    Parameters: cells, settings, results (as in `run_replications()`),
    replications (number of replications D; replications 1 to D are read).
    Returns: DataFrame with run, setting, N, method, gamma and R, R_sd, FDR,
    FDR_sd, TPR, TPR_sd; one row per Bayesian cell and three per EB cell.
    """
    def summary_row(cell, method, gamma, R, fdr, tpr):
        return {"run": cell["run"], "setting": cell["setting"], "N": cell["N"], "method": method,
                "gamma": gamma,
                "R": rmath.fround(rbase.r_mean_integer(R, na_rm=True), 2),
                "R_sd": rmath.fround(rbase.r_sd(R, na_rm=True), 2),
                "FDR": rmath.fround(rbase.r_mean(fdr, na_rm=True), 4),
                "FDR_sd": rmath.fround(rbase.r_sd(fdr, na_rm=True), 4),
                "TPR": rmath.fround(rbase.r_mean(tpr, na_rm=True), 4),
                "TPR_sd": rmath.fround(rbase.r_sd(tpr, na_rm=True), 4)}

    def zero_if_none(V, R):
        return np.where(R == 0, 0.0, _ratio(V, R))

    rows = []
    for cell in cells.to_dict("records"):
        files = [_replication_file(cell, d, settings, results) for d in range(1, replications + 1)]
        if not all(os.path.exists(file) for file in files):
            raise FileNotFoundError(f"replications missing for {os.path.dirname(files[0])}")
        reps = pd.concat([rbase.read_csv(file) for file in files], ignore_index=True)
        reps = {name: reps[name].to_numpy(dtype=np.float64) for name in reps.columns}
        if cell["method"] != "EB":
            rows.append(summary_row(cell, cell["method"], cell["gamma"], reps["R"],
                                    zero_if_none(reps["V"], reps["R"]),
                                    _ratio(reps["TP"], reps["n_alt"])))
            continue
        rows.append(summary_row(cell, "TS", np.nan, reps["TS_R"],
                                zero_if_none(reps["TS_V"], reps["TS_R"]),
                                _ratio(reps["TS_TP"], reps["n_alt"])))
        rows.append(summary_row(cell, "OS", np.nan, reps["OS_R"],
                                zero_if_none(reps["OS_V"], reps["OS_R"]),
                                _ratio(reps["OS_TP"], reps["n_alt"])))
        rows.append(summary_row(cell, "Storey", np.nan, reps["Storey_R"],
                                _ratio(reps["Storey_V"], reps["Storey_R"]),
                                _ratio(reps["Storey_TP"], reps["n_alt"])))
    return pd.DataFrame(rows)


def simulation_summaries(item, opts, settings, paths):
    """Summaries of the settings of a table or figure, by simulation.

    Runs the replications asked for that are not stored yet; when
    replications 1 to `opts["replications"]` are all stored, returns their
    summaries.

    Parameters: item ("table1", "table2", "figS1" or "figS3"), opts
    (options: `runs`, replications to run now; `replications`, D; `cores`),
    settings (settings of the paper), paths (project paths).
    Returns: DataFrame as `summarize_replications()`, or None if
    replications are still missing.
    """
    setup = item_settings(item)

    def with_method(method, gamma):
        return setup.assign(method=method, gamma=gamma)

    if item == "table2":
        cells = pd.concat([with_method(method, g) for g in PRIOR_GAMMAS.values()
                           for method in ("NP", "SP")], ignore_index=True)
    else:
        cells = pd.concat([with_method("NP", DEFAULT_GAMMA["NP"]),
                           with_method("GP", DEFAULT_GAMMA["GP"]),
                           with_method("SP", DEFAULT_GAMMA["SP"]),
                           with_method("EB", np.nan)], ignore_index=True)
    runs = parse_index_set(opts["runs"], count=True)
    run_replications(cells, [d for d in runs if d <= opts["replications"]], settings,
                     paths["results"], opts["cores"])
    done = all(os.path.exists(_replication_file(cell, d, settings, paths["results"]))
               for cell in cells.to_dict("records") for d in range(1, opts["replications"] + 1))
    if not done:
        log_message("replications 1 to %d are not all stored yet; no summary written",
                    opts["replications"])
        return None
    return summarize_replications(cells, opts["replications"], settings, paths["results"])


def _keys(table, columns):
    return list(zip(*(table[column].tolist() for column in columns)))


def archived_summaries(item, paths):
    """Archived summaries of the original runs for a table or figure.

    Parameters: item ("table1", "table2", "figS1" or "figS3"), paths
    (project paths).
    Returns: DataFrame with study, run, setting, N, method, gamma and the
    archived R, R_sd, FDR, FDR_sd, TPR, TPR_sd (missing where a run did not
    store them).
    """
    archive = rbase.read_csv(os.path.join(paths["archive"], "simulation_summaries.csv"))
    wanted = set(_keys(item_settings(item), ("run", "setting", "N")))
    archive = archive[[key in wanted for key in _keys(archive, ("run", "setting", "N"))]]
    study = "prior" if item == "table2" else "methods"
    return archive[archive["study"] == study].reset_index(drop=True)


def _match(keys, table_keys):
    """`match(keys, table_keys)`, zero-based, -1 where there is no match."""
    first = {}
    for position, key in enumerate(table_keys):
        first.setdefault(key, position)
    return [first.get(key, -1) for key in keys]


def _take(table, rows):
    """Rows of a table by position; a row of missing values for -1."""
    taken = table.iloc[[max(row, 0) for row in rows]].reset_index(drop=True)
    missing = [row < 0 for row in rows]
    if any(missing):
        taken = taken.astype({name: np.float64 for name in taken.columns
                              if pd.api.types.is_integer_dtype(taken[name])})
        taken.loc[missing, :] = np.nan
    return taken


def with_parameters(summaries):
    """Add the parameters of the settings to a table of summaries.

    Parameters: summaries (DataFrame with run, setting and N).
    Returns: `summaries` with the columns lambda, theta, pi0 and C.
    """
    every = simulation_settings()
    rows = _match(_keys(summaries, ("run", "setting", "N")), _keys(every, ("run", "setting", "N")))
    parameters = _take(every[["lambda", "theta", "pi0", "C"]], rows)
    return pd.concat([summaries.reset_index(drop=True), parameters], axis=1)


# ----- Layout of the tables -----

def summary_block(summaries, prefix):
    """The six summary columns of a block of a table, with a prefix.

    Parameters: summaries (rows of summaries, in the order of the table
    rows), prefix (prefix of the column names).
    Returns: DataFrame with <prefix>_R, _R_sd, _FDR, _FDR_sd, _TPR, _TPR_sd.
    """
    block = summaries[list(SUMMARY_COLUMNS)].reset_index(drop=True).astype(np.float64)
    return block.rename(columns={name: f"{prefix}_{name}" for name in SUMMARY_COLUMNS})


def table1_layout(summaries):
    """Table 1: R, FDR and TPR of the six methods by N, for C = 10 and C = 5.

    Parameters: summaries (summaries of the settings of
    `item_settings("table1")`).
    Returns: DataFrame with one row per N and method.
    """
    summaries = with_parameters(summaries)
    rows = [(N, method) for N in sorted(set(summaries["N"].tolist())) for method in METHODS]

    def block(C):
        found = _match([(N, method, C) for N, method in rows], _keys(summaries, ("N", "method", "C")))
        return _take(summaries, found)

    first = pd.DataFrame({"N": [N for N, _ in rows],
                          "method": [METHOD_LABELS[method] for _, method in rows]})
    return pd.concat([first, summary_block(block(10), "C10"), summary_block(block(5), "C5")], axis=1)


def table2_layout(summaries):
    """Table 2: R, FDR and TPR of NP and SP by N and prior.

    Parameters: summaries (summaries of the settings of
    `item_settings("table2")`).
    Returns: DataFrame with one row per N and prior.
    """
    rows = [(N, prior, gamma) for N in sorted(set(summaries["N"].tolist()))
            for prior, gamma in PRIOR_GAMMAS.items()]

    def block(method):
        found = _match([(N, gamma, method) for N, _, gamma in rows],
                       _keys(summaries, ("N", "gamma", "method")))
        return _take(summaries, found)

    first = pd.DataFrame({"N": [N for N, _, _ in rows], "prior": [prior for _, prior, _ in rows],
                          "gamma": [float(gamma) for _, _, gamma in rows]})
    return pd.concat([first, summary_block(block("NP"), "NP"), summary_block(block("SP"), "SP")],
                     axis=1)
