# -----------------------------------------------------------------------------
# File:    scripts/figS2_histograms.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Figure S2.
# Purpose: Histograms of nine simulated data sets with a ZIGP null (eta = 0.4,
#          theta = 0.3), a shifted binomial f1, pi0 = 0.8 and C = 5, for
#          N = 50, 100, 200 and lambda = 2, 3, 4: null counts in gray, non-null
#          counts in black.
# Inputs:  None (simulated data, seed 1439).
# Outputs: figures/FigS2<A-I>_simulated-counts_N<N>_lambda<lambda>.pdf,
#          results/FigS2_frequencies.csv (the plotted frequencies)
# Usage:   Rscript scripts/figS2_histograms.R
#          Always computed (--rerun changes nothing).
# Source:  barplots0504.R, rmixture() of simFDR_fcns.R.
# -----------------------------------------------------------------------------

source("R/setup.R")

opts <- parse_args(list())

frequencies <- do.call(rbind, lapply(seq_len(nrow(FIGS2_PANELS)), function(i) {
  setting <- data.frame(N = FIGS2_PANELS$N[i], eta = 0.4, lambda = FIGS2_PANELS$lambda[i],
                        theta = 0.3, pi0 = 0.8, f1 = "binomial", C = 5)
  data <- simulated_data(setting, PAPER$figS2_seed)
  plotted <- mixture_barplot(data$null, data$alt, output_file(PATHS$figures, figS2_file(i)),
                             PAPER$figure_size)
  data.frame(N = setting$N, lambda = setting$lambda, count = seq_along(plotted$null) - 1,
             null = plotted$null, non_null = plotted$alt)
}))
write.csv(frequencies, output_file(PATHS$results, "FigS2_frequencies.csv"), row.names = FALSE)
log_message("wrote %d figures to %s", nrow(FIGS2_PANELS), PATHS$figures)
