# -----------------------------------------------------------------------------
# File:    tests/testthat/test-plots.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Unit tests of R/plots.R: every plot function writes its PDF file
#          and returns the plotted values or the plot, whose layers hold the
#          expected marks; the names of the figure files.
# Inputs:  None (figures are written to temporary folders).
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/setup.R")
library(ggplot2)

#' Values as the script of Figure S1 passes them: two N for each method.
made_up_rates <- function() {
  data.frame(N = rep(c(50, 100), each = 6), method = rep(names(PROCEDURES), 2),
             value = seq(0.01, 0.12, by = 0.01))
}

test_that("frequency_barplot writes a pdf and returns the plotted frequencies", {
  file <- file.path(new_folder(), "bars.pdf")
  frequency <- count_frequencies(c(0, 0, 0, 1, 1, 4))
  plotted <- frequency_barplot(frequency, file, size = 3)
  expect_true(is_pdf(file))
  expect_equal(plotted, c("0" = 3, "1" = 2, "2" = 0, "3" = 0, "4" = 1))
  expect_null(dev.list())  # the device is closed
})

test_that("mixture_barplot writes a pdf and returns the two sets of frequencies", {
  file <- file.path(new_folder(), "mixture.pdf")
  plotted <- mixture_barplot(null = c(0, 0, 1, 2, 2, 2), alt = c(4, 6, 6), file, size = 3)
  expect_true(is_pdf(file))
  expect_equal(plotted$null, c("0" = 2, "1" = 1, "2" = 3, "3" = 0, "4" = 0, "5" = 0, "6" = 0))
  expect_equal(plotted$alt, c(0, 0, 0, 0, 1, 0, 2))
  expect_null(dev.list())
})

test_that("rate_plot draws a line and a point per value, and the level for the FDR", {
  folder <- new_folder()
  values <- made_up_rates()
  fdr <- rate_plot(values, "FDR", file.path(folder, "fdr.pdf"), size = 3)
  tpr <- rate_plot(values, "TPR", file.path(folder, "tpr.pdf"), size = 3)
  expect_true(is_pdf(file.path(folder, "fdr.pdf")))
  expect_true(is_pdf(file.path(folder, "tpr.pdf")))
  expect_s3_class(fdr, "ggplot")
  layers <- ggplot_build(fdr)$data
  expect_length(layers, 3)                       # lines, points, the line at 0.05
  expect_equal(nrow(layers[[1]]), 12)
  expect_equal(sort(layers[[2]]$y), values$value)
  expect_equal(sort(unique(layers[[2]]$x)), c(50, 100))
  expect_equal(length(unique(layers[[2]]$colour)), 6)
  expect_equal(layers[[3]]$yintercept, 0.05)
  expect_equal(fdr$theme$legend.position, "none")
  expect_length(ggplot_build(tpr)$data, 2)       # no line at 0.05, and the legend
  expect_false(identical(tpr$theme$legend.position, "none"))
  expect_equal(c(fdr$labels$y, tpr$labels$y), c("FDR", "TPR"))
  # the three Bayesian methods have solid lines, the others dashed
  lines <- unique(ggplot_build(tpr)$data[[1]][c("group", "linetype")])
  expect_equal(lines$linetype[order(lines$group)], rep(c("solid", "dashed"), each = 3))
})

test_that("fdr_table_plot draws one shaded cell per value", {
  file <- file.path(new_folder(), "table.pdf")
  fdr <- expand.grid(N = c(50, 100), lambda = c(1, 2), theta = c(0.15, 0.2, 0.3),
                     method = c("SP", "NP"), stringsAsFactors = FALSE)
  fdr$FDR <- seq(0.03, by = 0.002, length.out = nrow(fdr))
  plot <- fdr_table_plot(fdr, c("SP", "NP"), file)
  expect_true(is_pdf(file))
  expect_s3_class(plot, "ggplot")
  layers <- ggplot_build(plot)$data
  expect_equal(nrow(layers[[1]]), 24)
  expect_equal(sort(unique(as.numeric(layers[[1]]$x))), 1:6)   # method x theta
  expect_equal(sort(unique(as.numeric(layers[[1]]$y))), 1:4)   # N x lambda
  expect_setequal(layers[[2]]$label, sprintf("%.4f", fdr$FDR))
  # shades: below 0.05, 0.05 to 0.06, 0.06 to 0.07, 0.07 and above
  shade <- c("#FCE4D6", "#F4B084", "#C65911", "#833C0B")[findInterval(fdr$FDR, c(0.05, 0.06, 0.07)) + 1]
  expect_equal(as.vector(table(layers[[1]]$fill)[c("#FCE4D6", "#F4B084", "#C65911", "#833C0B")]),
               as.vector(table(factor(shade, c("#FCE4D6", "#F4B084", "#C65911", "#833C0B")))))
  expect_equal(levels(plot$data$column), c("SP 0.15", "SP 0.2", "SP 0.3", "NP 0.15", "NP 0.2", "NP 0.3"))
})

test_that("the figure files are named by panel", {
  expect_equal(fig1_file(1), "Fig1A_mutation-counts_cd00054.pdf")
  expect_equal(fig1_file(6), "Fig1F_mutation-counts_cd00063.pdf")
  expect_equal(figS1_file(1), "FigS1A_FDR_well-separated_setting16.pdf")
  expect_equal(figS1_file(4), "FigS1D_TPR_well-separated_setting18.pdf")
  expect_equal(figS2_file(1), "FigS2A_simulated-counts_N50_lambda2.pdf")
  expect_equal(figS2_file(9), "FigS2I_simulated-counts_N200_lambda4.pdf")
  expect_equal(nrow(FIGS1_PANELS), 4)
  expect_equal(nrow(FIGS2_PANELS), 9)
  expect_length(unique(vapply(1:9, figS2_file, character(1))), 9)
})
