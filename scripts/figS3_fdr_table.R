# -----------------------------------------------------------------------------
# File:    scripts/figS3_fdr_table.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Figure S3.
# Purpose: Table of the FDR of the six methods for N = 50, 100, 200, 500, 1000,
#          lambda = 1, 2, 3, 4 and theta = 0.15, 0.2, 0.3 (ZIGP null with
#          eta = 0.4, shifted binomial f1, pi0 = 0.8, C = 5), shaded by FDR.
# Inputs:  archive/simulation_summaries.csv (FDR of the original runs run0219,
#          run0227 and run0207); none with --rerun (simulated data).
# Outputs: results/FigS3_fdr_table.csv (rows N and lambda; columns method and
#          theta, in the order of the published table),
#          figures/FigS3_fdr_table.pdf; with --rerun also
#          results/replications/<setting>/<method>/repDDDD.csv.
# Usage:   Rscript scripts/figS3_fdr_table.R [--rerun] [--runs=1000]
#          [--replications=1000] [--cores=4]
#          Options as in scripts/table1_simulation.R.
# Source:  facet() of simFDR_fcns.R: sim_initial0219.R settings 1 to 30,
#          sim_initial0227.R settings 1 to 20, initial_vals0207.R settings 10
#          and 14; sheet "Full Format" of Heatmaps AOAS.xlsx (layout).
# -----------------------------------------------------------------------------

source("R/setup.R")

opts <- parse_args(simulation_options(PAPER))
TABLE_METHODS <- c("SP", "GP", "NP", "OS", "Storey", "TS")  # column order of the published table

summaries <- if (opts$rerun) {
  simulation_summaries("figS3", opts, PAPER, PATHS)
} else {
  archived_summaries("figS3", PATHS)
}
if (!is.null(summaries)) {
  fdr <- with_parameters(summaries)[c("N", "lambda", "theta", "method", "FDR")]
  rows <- unique(fdr[order(fdr$N, fdr$lambda), c("N", "lambda")])
  columns <- expand.grid(theta = sort(unique(fdr$theta)), method = TABLE_METHODS,
                         stringsAsFactors = FALSE)
  table <- cbind(rows, sapply(seq_len(nrow(columns)), function(k) {
    fdr$FDR[match(paste(rows$N, rows$lambda, columns$theta[k], columns$method[k]),
                  paste(fdr$N, fdr$lambda, fdr$theta, fdr$method))]
  }))
  names(table) <- c("N", "lambda", sprintf("%s_theta%s", columns$method, columns$theta))
  name <- study_name("FigS3_fdr_table", opts, PAPER)
  file <- output_file(PATHS$results, paste0(name, ".csv"))
  write.csv(table, file, row.names = FALSE)
  fdr_table_plot(fdr, TABLE_METHODS, output_file(PATHS$figures, paste0(name, ".pdf")))
  log_message("wrote %s", file)
}
