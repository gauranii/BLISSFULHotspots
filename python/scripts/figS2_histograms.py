# -----------------------------------------------------------------------------
# File:    python/scripts/figS2_histograms.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Figure S2.
# Purpose: Bar plots of nine simulated data sets (N = 50, 100, 200; lambda = 2, 3, 4), null and non-null counts.
# Inputs:  None (simulated data, seed 1439).
# Outputs: figures/FigS2<A-I>_simulated-counts_N<N>_lambda<lambda>.pdf, results/FigS2_frequencies.csv.
# Usage:   python python/scripts/figS2_histograms.py 
#          Options and outputs as scripts/figS2_histograms.R, which this file follows.
# Source:  scripts/figS2_histograms.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from localfdr.plots import FIGS2_PANELS, figS2_file, mixture_barplot  # noqa: E402
from localfdr.rbase import write_csv  # noqa: E402
from localfdr.setup import PAPER, PATHS, log_message, output_file, parse_args  # noqa: E402
from localfdr.simulate import simulated_data  # noqa: E402

opts = parse_args({})

parts = []
for i, (N, lam) in enumerate(zip(FIGS2_PANELS["N"], FIGS2_PANELS["lambda"]), start=1):
    setting = {"N": N, "eta": 0.4, "lambda": lam, "theta": 0.3, "pi0": 0.8, "f1": "binomial", "C": 5}
    data = simulated_data(setting, PAPER["figS2_seed"])
    plotted = mixture_barplot(data["null"], data["alt"],
                              output_file(PATHS["figures"], figS2_file(i)), PAPER["figure_size"])
    parts.append(pd.DataFrame({"N": N, "lambda": lam, "count": np.arange(len(plotted["null"])),
                               "null": plotted["null"], "non_null": plotted["alt"]}))
write_csv(pd.concat(parts, ignore_index=True), output_file(PATHS["results"], "FigS2_frequencies.csv"))
log_message("wrote %d figures to %s", len(FIGS2_PANELS["N"]), PATHS["figures"])
