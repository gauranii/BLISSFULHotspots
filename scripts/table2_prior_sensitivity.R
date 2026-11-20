# -----------------------------------------------------------------------------
# File:    scripts/table2_prior_sensitivity.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Table 2.
# Purpose: Prior sensitivity of the nonparametric and semiparametric methods:
#          R, FDR and TPR (with standard deviations) for five Dirichlet
#          concentrations gamma = 0.02 (Objective), 0.5 (Jeffreys'), 1
#          (Uniform), 1.5 (Near Uniform) and 3 (Informative), in the settings
#          of Table 1 with C = 10.
# Inputs:  archive/simulation_summaries.csv (summaries of the original run);
#          none with --rerun (simulated data).
# Outputs: results/Table2_prior-sensitivity.csv (the table); with --rerun also
#          results/replications/<setting>/<method>/repDDDD.csv.
# Usage:   Rscript scripts/table2_prior_sensitivity.R [--rerun] [--runs=1000]
#          [--replications=1000] [--cores=4]
#          Options as in scripts/table1_simulation.R. The replications with
#          gamma = 1 are those of Table 1 and are shared with it.
# Source:  prior(hyper_vec, j) of simFDR_fcns.R (run0224 version), settings
#          1, 2, 3, 6, 11 of sim_initial0225.R.
# -----------------------------------------------------------------------------

source("R/setup.R")

opts <- parse_args(simulation_options(PAPER))

summaries <- if (opts$rerun) {
  simulation_summaries("table2", opts, PAPER, PATHS)
} else {
  archived_summaries("table2", PATHS)
}
if (!is.null(summaries)) {
  file <- output_file(PATHS$results, paste0(study_name("Table2_prior-sensitivity", opts, PAPER), ".csv"))
  write.csv(table2_layout(summaries), file, row.names = FALSE)
  log_message("wrote %s", file)
}
