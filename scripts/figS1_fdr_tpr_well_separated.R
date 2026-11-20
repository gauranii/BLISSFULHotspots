# -----------------------------------------------------------------------------
# File:    scripts/figS1_fdr_tpr_well_separated.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Figure S1.
# Purpose: FDR and TPR of the six methods against N = 50, 100, ..., 1000 when
#          f0 and f1 are well separated: zero-inflated Poisson null (eta = 0.4,
#          lambda = 1.5), shifted binomial f1, C = 10, with pi0 = 0.35 (panels
#          A, B) and pi0 = 0.80 (panels C, D).
# Inputs:  archive/simulation_summaries.csv (summaries of the original run
#          run0207, settings 16 and 18); none with --rerun (simulated data).
# Outputs: figures/FigS1<A-D>_<FDR|TPR>_well-separated_setting<16|18>.pdf,
#          results/FigS1_values.csv (the plotted values); with --rerun also
#          results/replications/<setting>/<method>/repDDDD.csv.
# Usage:   Rscript scripts/figS1_fdr_tpr_well_separated.R [--rerun] [--runs=1000]
#          [--replications=1000] [--cores=4]
#          Options as in scripts/table1_simulation.R.
# Source:  facet(n_facet, j) of simFDR_fcns.R (run0207 version), settings 16
#          and 18 of initial_vals0207.R; facet_fdr() and facet_tpr() of the
#          plotting run of 2021-02-17.
# -----------------------------------------------------------------------------

source("R/setup.R")

opts <- parse_args(simulation_options(PAPER))

summaries <- if (opts$rerun) {
  simulation_summaries("figS1", opts, PAPER, PATHS)
} else {
  archived_summaries("figS1", PATHS)
}
if (!is.null(summaries)) {
  values <- do.call(rbind, lapply(seq_len(nrow(FIGS1_PANELS)), function(i) {
    panel <- FIGS1_PANELS[i, ]
    rows <- summaries[summaries$setting == panel$setting, ]
    plotted <- data.frame(setting = panel$setting, measure = panel$measure, N = rows$N,
                          method = rows$method, value = rows[[panel$measure]])
    name <- sub("FigS1", study_name("FigS1", opts, PAPER), figS1_file(i), fixed = TRUE)
    rate_plot(plotted, panel$measure, output_file(PATHS$figures, name), PAPER$figure_size)
    plotted
  }))
  file <- output_file(PATHS$results, paste0(study_name("FigS1_values", opts, PAPER), ".csv"))
  write.csv(values, file, row.names = FALSE)
  log_message("wrote %s and %d figures", file, nrow(FIGS1_PANELS))
}
