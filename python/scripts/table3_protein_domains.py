# -----------------------------------------------------------------------------
# File:    python/scripts/table3_protein_domains.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Table 3.
# Purpose: Number of positions declared hotspots in the protein domains, by method.
# Inputs:  data/domains/<domain>.txt.
# Outputs: results/Table3_protein-domains.csv (with --domains: Table3_protein-domains_<domains>.csv).
# Usage:   python python/scripts/table3_protein_domains.py [--cores=4] [--domains=all] [--poisson]
#          Options and outputs as scripts/table3_protein_domains.R, which this file follows.
# Source:  scripts/table3_protein_domains.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from localfdr.rbase import write_csv  # noqa: E402
from localfdr.real_data import DOMAINS, domain_rejections  # noqa: E402
from localfdr.setup import PAPER, PATHS, log_message, output_file, parse_args  # noqa: E402

opts = parse_args({"cores": 4, "domains": "all", "poisson": False})
domains = sorted(DOMAINS) if opts["domains"] == "all" else opts["domains"].split(",")
models = ["NP", "SP"] + (["P"] if opts["poisson"] else []) + ["GP"]

table = domain_rejections(domains, models, PAPER, PATHS, opts["cores"])
name = "Table3_protein-domains.csv" if opts["domains"] == "all" else \
    "Table3_protein-domains_%s.csv" % "-".join(domains)
file = output_file(PATHS["results"], name)
write_csv(table, file)
print(table.to_string())
log_message("wrote %s", file)
