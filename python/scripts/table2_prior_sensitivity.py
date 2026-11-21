# -----------------------------------------------------------------------------
# File:    python/scripts/table2_prior_sensitivity.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Table 2.
# Purpose: R, FDR and TPR of the NP and SP methods by N and Dirichlet prior.
# Inputs:  archive/simulation_summaries.csv; none with --rerun (simulated data).
# Outputs: results/Table2_prior-sensitivity.csv; with --rerun also results/replications/.
# Usage:   python python/scripts/table2_prior_sensitivity.py [--rerun] [--runs=1000] [--replications=1000] [--cores=4]
#          Options and outputs as scripts/table2_prior_sensitivity.R, which this file follows.
# Source:  scripts/table2_prior_sensitivity.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from localfdr.rbase import write_csv  # noqa: E402
from localfdr.setup import (PAPER, PATHS, log_message, output_file, parse_args,  # noqa: E402
                            simulation_options, study_name)
from localfdr.simulate import archived_summaries, simulation_summaries, table2_layout  # noqa: E402

opts = parse_args(simulation_options(PAPER))

if opts["rerun"]:
    summaries = simulation_summaries("table2", opts, PAPER, PATHS)
else:
    summaries = archived_summaries("table2", PATHS)
if summaries is not None:
    file = output_file(PATHS["results"],
                       study_name("Table2_prior-sensitivity", opts, PAPER) + ".csv")
    write_csv(table2_layout(summaries), file)
    log_message("wrote %s", file)
