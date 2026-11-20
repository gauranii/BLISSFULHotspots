# -----------------------------------------------------------------------------
# File:    scripts/table1_simulation.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Table 1.
# Purpose: Average number of rejections R, FDR and TPR (with standard
#          deviations) of six methods in simulated data with a ZIGP null
#          (eta = 0.4, lambda = 2, theta = 0.3), a shifted binomial f1 and
#          pi0 = 0.8, for N = 50, 100, 200, 500, 1000 and C = 10 and C = 5.
# Inputs:  archive/simulation_summaries.csv (summaries of the original runs
#          run0224, C = 10, and run0219, C = 5); none with --rerun (simulated
#          data).
# Outputs: results/Table1_simulation.csv (the table); with --rerun also
#          results/replications/<setting>/<method>/repDDDD.csv.
# Usage:   Rscript scripts/table1_simulation.R [--rerun] [--runs=1000]
#          [--replications=1000] [--cores=4]
#          --rerun: run the replications here instead of reading the archive.
#          --runs: N (replications 1 to N), a:b or a comma list; finished
#            replications are kept.
#          --replications: replications D of the study; the table is written
#            when replications 1 to D are all stored (with D other than 1000
#            as Table1_simulation_first<D>-replications.csv).
# Source:  run0224: sim_initial0225.R settings 1, 2, 3, 6, 11; run0219:
#          sim_initial0219.R settings 23, 5, 29, 11, 17; facet() of
#          simFDR_fcns.R.
# -----------------------------------------------------------------------------

source("R/setup.R")

opts <- parse_args(simulation_options(PAPER))

summaries <- if (opts$rerun) {
  simulation_summaries("table1", opts, PAPER, PATHS)
} else {
  archived_summaries("table1", PATHS)
}
if (!is.null(summaries)) {
  file <- output_file(PATHS$results, paste0(study_name("Table1_simulation", opts, PAPER), ".csv"))
  write.csv(table1_layout(summaries), file, row.names = FALSE)
  log_message("wrote %s", file)
}
