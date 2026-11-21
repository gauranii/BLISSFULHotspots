# -----------------------------------------------------------------------------
# File:    python/scripts/figS1_fdr_tpr_well_separated.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Figure S1.
# Purpose: FDR and TPR against N for the six methods in two well-separated settings.
# Inputs:  archive/simulation_summaries.csv; none with --rerun (simulated data).
# Outputs: figures/FigS1<A-D>_<measure>_well-separated_setting<k>.pdf, results/FigS1_values.csv.
# Usage:   python python/scripts/figS1_fdr_tpr_well_separated.py [--rerun] [--runs=1000] [--replications=1000] [--cores=4]
#          Options and outputs as scripts/figS1_fdr_tpr_well_separated.R, which this file follows.
# Source:  scripts/figS1_fdr_tpr_well_separated.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import pandas as pd  # noqa: E402

from localfdr.plots import FIGS1_PANELS, figS1_file, rate_plot  # noqa: E402
from localfdr.rbase import write_csv  # noqa: E402
from localfdr.setup import (PAPER, PATHS, log_message, output_file, parse_args,  # noqa: E402
                            simulation_options, study_name)
from localfdr.simulate import archived_summaries, simulation_summaries  # noqa: E402

opts = parse_args(simulation_options(PAPER))

if opts["rerun"]:
    summaries = simulation_summaries("figS1", opts, PAPER, PATHS)
else:
    summaries = archived_summaries("figS1", PATHS)
if summaries is not None:
    parts = []
    for i, (setting, measure) in enumerate(zip(FIGS1_PANELS["setting"], FIGS1_PANELS["measure"]),
                                           start=1):
        rows = summaries[summaries["setting"] == setting]
        plotted = pd.DataFrame({"setting": setting, "measure": measure, "N": rows["N"].to_numpy(),
                                "method": rows["method"].to_numpy(),
                                "value": rows[measure].to_numpy(dtype=float)})
        name = figS1_file(i).replace("FigS1", study_name("FigS1", opts, PAPER), 1)
        rate_plot(plotted, measure, output_file(PATHS["figures"], name), PAPER["figure_size"])
        parts.append(plotted)
    file = output_file(PATHS["results"], study_name("FigS1_values", opts, PAPER) + ".csv")
    write_csv(pd.concat(parts, ignore_index=True), file)
    log_message("wrote %s and %d figures", file, len(FIGS1_PANELS["setting"]))
