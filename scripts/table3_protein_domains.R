# -----------------------------------------------------------------------------
# File:    scripts/table3_protein_domains.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Table 3.
# Purpose: Number of positions declared mutational hotspots in six protein
#          domains by the empirical Bayes method (EB) and by the Bayesian
#          methods with a nonparametric (NP), semiparametric (SP) and
#          generalized Poisson (GP) alternative, at local fdr 0.05.
# Inputs:  data/domains/<domain>.txt
# Outputs: results/Table3_protein-domains.csv
# Usage:   Rscript scripts/table3_protein_domains.R [--cores=4] [--domains=all]
#          [--poisson]
#          --domains: comma list of domain names.
#          --poisson: also run the Poisson alternative (column P).
#          Always computed from the data (--rerun changes nothing): one chain
#          of 11000 iterations per domain and method, seed 1234.
# Source:  real_scenarios_wz.R, real_scenarios.R (gibbs_out, RealDataFDR).
# -----------------------------------------------------------------------------

source("R/setup.R")

opts <- parse_args(list(cores = 4, domains = "all", poisson = FALSE))
domains <- if (opts$domains == "all") sort(DOMAINS) else strsplit(opts$domains, ",")[[1]]
models <- c("NP", "SP", if (opts$poisson) "P", "GP")

table <- domain_rejections(domains, models, PAPER, PATHS, opts$cores)
name <- if (opts$domains == "all") "Table3_protein-domains.csv" else
  sprintf("Table3_protein-domains_%s.csv", paste(domains, collapse = "-"))
file <- output_file(PATHS$results, name)
write.csv(table, file, row.names = FALSE)
print(table)
log_message("wrote %s", file)
