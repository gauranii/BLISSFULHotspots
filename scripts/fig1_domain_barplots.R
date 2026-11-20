# -----------------------------------------------------------------------------
# File:    scripts/fig1_domain_barplots.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Figure 1.
# Purpose: Bar plots of six protein domains: the number of positions at each
#          mutation count.
# Inputs:  data/domains/<domain>.txt
# Outputs: figures/Fig1<A-F>_mutation-counts_<domain>.pdf,
#          results/Fig1_frequencies.csv (the plotted frequencies)
# Usage:   Rscript scripts/fig1_domain_barplots.R
#          Always computed from the data (--rerun changes nothing). As in the
#          published figure, the bars include position 0 of each file.
# Source:  protein.R.
# -----------------------------------------------------------------------------

source("R/setup.R")

opts <- parse_args(list())

frequencies <- do.call(rbind, lapply(seq_along(DOMAINS), function(i) {
  frequency <- count_frequencies(read_domain(DOMAINS[i], PATHS$domains)$with_position_0)
  frequency_barplot(frequency, output_file(PATHS$figures, fig1_file(i)), PAPER$figure_size)
  data.frame(domain = DOMAINS[i], count = as.integer(names(frequency)), positions = frequency)
}))
write.csv(frequencies, output_file(PATHS$results, "Fig1_frequencies.csv"), row.names = FALSE)
log_message("wrote %d figures to %s", length(DOMAINS), PATHS$figures)
