# -----------------------------------------------------------------------------
# File:    python/scripts/fig1_domain_barplots.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Figure 1.
# Purpose: Bar plots of the number of positions at each mutation count, one per domain.
# Inputs:  data/domains/<domain>.txt.
# Outputs: figures/Fig1<A-F>_mutation-counts_<domain>.pdf, results/Fig1_frequencies.csv.
# Usage:   python python/scripts/fig1_domain_barplots.py 
#          Options and outputs as scripts/fig1_domain_barplots.R, which this file follows.
# Source:  scripts/fig1_domain_barplots.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import pandas as pd  # noqa: E402

from localfdr.plots import fig1_file, frequency_barplot  # noqa: E402
from localfdr.rbase import write_csv  # noqa: E402
from localfdr.real_data import DOMAINS, count_frequencies, read_domain  # noqa: E402
from localfdr.setup import PAPER, PATHS, log_message, output_file, parse_args  # noqa: E402

opts = parse_args({})

parts = []
for i, domain in enumerate(DOMAINS, start=1):
    frequency = count_frequencies(read_domain(domain, PATHS["domains"])["with_position_0"])
    frequency_barplot(frequency, output_file(PATHS["figures"], fig1_file(i)), PAPER["figure_size"])
    parts.append(pd.DataFrame({"domain": domain, "count": list(frequency),
                               "positions": list(frequency.values())}))
write_csv(pd.concat(parts, ignore_index=True), output_file(PATHS["results"], "Fig1_frequencies.csv"))
log_message("wrote %d figures to %s", len(DOMAINS), PATHS["figures"])
