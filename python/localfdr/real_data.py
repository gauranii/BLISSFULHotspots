# -----------------------------------------------------------------------------
# File:    python/localfdr/real_data.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Table 3, Figure 1.
# Purpose: Reader of the protein domain files and the analysis of a domain:
#          the number of positions declared mutational hotspots by the
#          empirical Bayes method and by the Bayesian methods.
# Inputs:  None.
# Outputs: None.
# Usage:   from localfdr.real_data import read_domain, domain_rejections
# Source:  R/real_data.R.
# -----------------------------------------------------------------------------
import os

import numpy as np
import pandas as pd

from . import rbase, rng
from .bayes_fdr import bayes_local_fdr, sampler_settings
from .eb_fdr import eb_rejections

# The six protein domains, in the order of Figure 1.
DOMAINS = ("cd00054", "cd00083", "cd00590", "smart00179", "smart00209", "cd00063")


def read_domain(domain, folder):
    """Read the mutation counts of a protein domain.

    A domain file has the number of positions N on its first line, then one
    line `position:count:mutations` per position, starting with a position 0
    that holds no mutation.

    Parameters: domain (name of the domain, for example "cd00083"), folder
    (folder of the domain files).
    Returns: dict with `N`, `counts` (mutation counts of positions 1 to N,
    integers as in R) and `with_position_0` (the same with the count of
    position 0 in front).
    """
    with open(os.path.join(folder, domain + ".txt")) as f:
        lines = [line.rstrip("\n") for line in f if line.strip()]
    counts = np.array([int(line.split(":")[1]) for line in lines[1:]], dtype=np.int64)
    return {"N": int(lines[0].split(":")[0]), "counts": counts[1:], "with_position_0": counts}


def count_frequencies(x):
    """Number of positions at each mutation count 0, 1, ..., max.

    Parameters: x (mutation counts).
    Returns: dict from count to frequency, in increasing order of the count.
    """
    x = np.asarray(x, dtype=np.float64)
    largest = int(x.max())
    frequency = rbase.count_frequencies_from_zero(x, largest)
    return {count: int(frequency[count]) for count in range(largest + 1)}


def domain_bayes_rejections(x, model, settings):
    """Hotspot positions of a domain by one Bayesian method.

    Runs the sampler as in the analysis of the paper (seed, iterations and
    Dirichlet concentration from `settings`) and counts the positions whose
    mutation count has posterior mean local fdr at most alpha.

    Parameters: x (mutation counts of the positions), model ("NP", "SP", "P"
    or "GP"), settings (settings of the paper).
    Returns: number of positions declared hotspots (float).
    """
    rng.set_seed(settings["real_seed"])
    fit = bayes_local_fdr(x, model, settings["real_gamma"][model],
                          sampler_settings(lambda_init=1, theta_init=0.5),
                          draws=settings["real_draws"], burn=settings["real_burn"],
                          version="2018-03-25")
    return rbase.r_sum(fit["frequency"][fit["fdr"] <= settings["alpha"]])


def domain_rejections(domains, models, settings, paths, cores=1):
    """Hotspot positions of the domains by every method (Table 3).

    Parameters: domains (names of the domains), models (Bayesian models to
    run), settings (settings of the paper), paths (project paths), cores
    (number of parallel processes).
    Returns: DataFrame with domain, N, EB and one column per model.
    """
    domains, models = list(domains), list(models)
    jobs = [(model, domain) for domain in domains for model in models]

    def run(job):
        model, domain = job
        x = read_domain(domain, paths["domains"])["counts"]
        return domain_bayes_rejections(x, model, settings)

    counts = rbase.parallel_map(run, jobs, cores)
    failed = [count for count in counts if not isinstance(count, float)]
    if failed:
        raise RuntimeError(f"a chain failed: {failed[0]!r}")
    bayes = np.array(counts).reshape(len(domains), len(models))
    table = {"domain": domains,
             "N": [float(len(read_domain(d, paths["domains"])["counts"])) for d in domains],
             "EB": [eb_rejections(read_domain(d, paths["domains"])["counts"], settings["alpha"])
                    for d in domains]}
    for j, model in enumerate(models):
        table[model] = bayes[:, j]
    return pd.DataFrame(table)
