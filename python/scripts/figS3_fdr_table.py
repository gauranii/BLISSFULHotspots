# -----------------------------------------------------------------------------
# File:    python/scripts/figS3_fdr_table.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Figure S3.
# Purpose: Shaded table of the FDR of the six methods by N, lambda and theta.
# Inputs:  archive/simulation_summaries.csv; none with --rerun (simulated data).
# Outputs: results/FigS3_fdr_table.csv, figures/FigS3_fdr_table.pdf.
# Usage:   python python/scripts/figS3_fdr_table.py [--rerun] [--runs=1000] [--replications=1000] [--cores=4]
#          Options and outputs as scripts/figS3_fdr_table.R, which this file follows.
# Source:  scripts/figS3_fdr_table.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import pandas as pd  # noqa: E402

from localfdr.plots import fdr_table_plot  # noqa: E402
from localfdr.rbase import format_number, write_csv  # noqa: E402
from localfdr.setup import (PAPER, PATHS, log_message, output_file, parse_args,  # noqa: E402
                            simulation_options, study_name)
from localfdr.simulate import archived_summaries, simulation_summaries, with_parameters  # noqa: E402

opts = parse_args(simulation_options(PAPER))
TABLE_METHODS = ["SP", "GP", "NP", "OS", "Storey", "TS"]  # column order of the published table

if opts["rerun"]:
    summaries = simulation_summaries("figS3", opts, PAPER, PATHS)
else:
    summaries = archived_summaries("figS3", PATHS)
if summaries is not None:
    fdr = with_parameters(summaries)[["N", "lambda", "theta", "method", "FDR"]]
    rows = sorted(set(zip(fdr["N"].tolist(), fdr["lambda"].tolist())))
    first = {}  # the first FDR of every (N, lambda, theta, method), as match() finds it
    for N, lam, theta, method, value in zip(fdr["N"], fdr["lambda"], fdr["theta"], fdr["method"],
                                            fdr["FDR"]):
        first.setdefault((N, lam, theta, method), float(value))
    table = {"N": [N for N, _ in rows], "lambda": [lam for _, lam in rows]}
    for method in TABLE_METHODS:
        for theta in sorted(set(fdr["theta"].tolist())):
            table["%s_theta%s" % (method, format_number(theta))] = [
                first.get((N, lam, theta, method), float("nan")) for N, lam in rows]
    name = study_name("FigS3_fdr_table", opts, PAPER)
    file = output_file(PATHS["results"], name + ".csv")
    write_csv(pd.DataFrame(table), file)
    fdr_table_plot(fdr, TABLE_METHODS, output_file(PATHS["figures"], name + ".pdf"))
    log_message("wrote %s", file)
