# -----------------------------------------------------------------------------
# File:    python/scripts/table1_simulation.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Table 1.
# Purpose: R, FDR and TPR of the six methods by N, for C = 10 and C = 5.
# Inputs:  archive/simulation_summaries.csv; none with --rerun (simulated data).
# Outputs: results/Table1_simulation.csv; with --rerun also results/replications/.
# Usage:   python python/scripts/table1_simulation.py [--rerun] [--runs=1000] [--replications=1000] [--cores=4]
#          Options and outputs as scripts/table1_simulation.R, which this file follows.
# Source:  scripts/table1_simulation.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from localfdr.rbase import write_csv  # noqa: E402
from localfdr.setup import (PAPER, PATHS, log_message, output_file, parse_args,  # noqa: E402
                            simulation_options, study_name)
from localfdr.simulate import archived_summaries, simulation_summaries, table1_layout  # noqa: E402

opts = parse_args(simulation_options(PAPER))

if opts["rerun"]:
    summaries = simulation_summaries("table1", opts, PAPER, PATHS)
else:
    summaries = archived_summaries("table1", PATHS)
if summaries is not None:
    file = output_file(PATHS["results"], study_name("Table1_simulation", opts, PAPER) + ".csv")
    write_csv(table1_layout(summaries), file)
    log_message("wrote %s", file)
