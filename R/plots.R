# -----------------------------------------------------------------------------
# File:    R/plots.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Figures 1, S1, S2 and S3.
# Purpose: The plots of the paper: bar plots of count frequencies, line plots
#          of FDR and TPR against N, and the shaded table of FDR values.
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R")
# Source:  protein.R (Figure 1), barplots0504.R (Figure S2), facet_fdr and
#          facet_tpr of simFDR_fcns.R (Figure S1), the sheet "Full Format" of
#          Heatmaps AOAS.xlsx (Figure S3).
# -----------------------------------------------------------------------------

# Names of the six methods in the legend of Figure S1 (drawn in this order).
PROCEDURES <- c(NP = "Bayesian - NP", GP = "Bayesian - P", SP = "Bayesian - SP",
                OS = "EB - One Stage", TS = "EB - Two Stage", Storey = "Storey")

#' Bar plot of the number of positions at each count (Figure 1).
#'
#' @param frequency Named vector of frequencies (`count_frequencies()`).
#' @param file Path of the PDF to write.
#' @param size Width and height in inches.
#' @return The plotted frequencies, invisibly.
frequency_barplot <- function(frequency, file, size) {
  pdf(file, width = size, height = size)
  on.exit(dev.off())
  barplot(frequency)
  invisible(frequency)
}

#' Bar plot of a simulated data set: null counts in gray, non-null counts in
#' black, and a dashed line right above the cutoff C = 5 (Figure S2).
#'
#' @param null,alt Null and non-null counts.
#' @param file Path of the PDF to write.
#' @param size Width and height in inches.
#' @return List with the plotted frequencies `null` and `alt`, invisibly.
mixture_barplot <- function(null, alt, file, size) {
  K <- max(alt)
  null_frequency <- c(sum(null == 0), tabulate(null, K))
  alt_frequency <- c(sum(alt == 0), tabulate(alt, K))
  names(null_frequency) <- c(0, seq(1, K, 1))
  pdf(file, width = size, height = size)
  on.exit(dev.off())
  barplot(null_frequency, cex.axis = 0.1, col = "gray")
  barplot(alt_frequency, col = "black", add = TRUE)
  abline(v = 7.2, lwd = 1.5, lty = 2, col = "blue")
  invisible(list(null = null_frequency, alt = alt_frequency))
}

#' Line plot of FDR or TPR against N for the six methods (Figure S1).
#'
#' @param values Data frame with N, method (names of `PROCEDURES`) and value.
#' @param measure "FDR" (with a line at 0.05 and no legend) or "TPR" (with
#'   the legend).
#' @param file Path of the PDF to write.
#' @param size Width and height in inches.
#' @return The plot (a ggplot object), invisibly.
rate_plot <- function(values, measure, file, size) {
  library(ggplot2)
  values$Procedure <- PROCEDURES[values$method]
  plot <- ggplot(values, aes(x = N, y = value, group = Procedure)) +
    geom_line(aes(linetype = Procedure, colour = Procedure), linewidth = 0.35) +
    scale_linetype_manual(values = c("solid", "solid", "solid", "dashed", "dashed", "dashed")) +
    xlab("N") + ylab(measure) +
    geom_point(aes(shape = Procedure, colour = Procedure), size = 1) +
    scale_shape_manual(values = c(8, 15, 17, 1, 2, 3)) +
    scale_color_brewer(palette = "Dark2")
  if (measure == "FDR") {
    plot <- plot + geom_hline(yintercept = 0.05, linewidth = 0.5, linetype = "solid") +
      theme_bw() + theme(legend.position = "none")
  } else {
    plot <- plot + theme_bw()
  }
  pdf(file, width = size, height = size)
  on.exit(dev.off())
  print(plot)
  invisible(plot)
}

#' Shaded table of FDR values by N, lambda, method and theta (Figure S3).
#'
#' Cells are shaded by FDR: below 0.05, 0.05 to 0.06, 0.06 to 0.07, 0.07 and
#' above.
#'
#' @param fdr Data frame with N, lambda, theta, method and FDR.
#' @param methods Methods in the order of the columns.
#' @param file Path of the PDF to write.
#' @return The plot (a ggplot object), invisibly.
fdr_table_plot <- function(fdr, methods, file) {
  library(ggplot2)
  fdr$column <- factor(paste(fdr$method, fdr$theta), levels = paste(
    rep(methods, each = 3), sort(unique(fdr$theta))))
  fdr$row <- factor(paste(fdr$N, fdr$lambda), levels = rev(unique(paste(fdr$N, fdr$lambda)[
    order(fdr$N, fdr$lambda)])))
  fdr$shade <- cut(fdr$FDR, c(-Inf, 0.05, 0.06, 0.07, Inf), right = FALSE,
                   labels = c("below 0.05", "0.05 to 0.06", "0.06 to 0.07", "0.07 and above"))
  plot <- ggplot(fdr, aes(x = column, y = row, fill = shade)) +
    geom_tile(colour = "black") +
    geom_text(aes(label = sprintf("%.4f", FDR), colour = shade), size = 2.6, show.legend = FALSE) +
    scale_fill_manual(values = c("#FCE4D6", "#F4B084", "#C65911", "#833C0B"), drop = FALSE,
                      name = "FDR") +
    scale_colour_manual(values = c("black", "black", "black", "white"), drop = FALSE) +
    scale_x_discrete(position = "top", labels = function(l) sub(" ", "\ntheta = ", l)) +
    scale_y_discrete(labels = function(l) sub("^(\\S+) (\\S+)$", "N = \\1, lambda = \\2", l)) +
    xlab(NULL) + ylab(NULL) + theme_minimal() +
    theme(panel.grid = element_blank(), legend.position = "bottom")
  pdf(file, width = 14, height = 6.5)
  on.exit(dev.off())
  print(plot)
  invisible(plot)
}

# ----- Names of the figure files -----

# Panels of Figure S1: setting of run0207 and measure.
FIGS1_PANELS <- data.frame(setting = c(16, 16, 18, 18), measure = c("FDR", "TPR", "FDR", "TPR"))
# Panels of Figure S2: N and lambda.
FIGS2_PANELS <- data.frame(N = rep(c(50, 100, 200), 3), lambda = rep(c(2, 3, 4), each = 3))

#' File name of a panel of Figure 1.
#' @param i Panel number (domain `DOMAINS[i]`).
#' @return File name.
fig1_file <- function(i) sprintf("Fig1%s_mutation-counts_%s.pdf", LETTERS[i], DOMAINS[i])

#' File name of a panel of Figure S1.
#' @param i Panel number (row of `FIGS1_PANELS`).
#' @return File name.
figS1_file <- function(i) {
  sprintf("FigS1%s_%s_well-separated_setting%d.pdf", LETTERS[i], FIGS1_PANELS$measure[i],
          FIGS1_PANELS$setting[i])
}

#' File name of a panel of Figure S2.
#' @param i Panel number (row of `FIGS2_PANELS`).
#' @return File name.
figS2_file <- function(i) {
  sprintf("FigS2%s_simulated-counts_N%d_lambda%d.pdf", LETTERS[i], FIGS2_PANELS$N[i],
          FIGS2_PANELS$lambda[i])
}
