# -----------------------------------------------------------------------------
# File:    archive/reduce_archive.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Tables 1 and 2, Figures S1 and S3.
# Purpose: Copy from the output folders of the original simulation runs the
#          values the table and figure scripts read, into the csv files kept
#          in archive/, and copy the published figure files under the names
#          of their figures.
# Inputs:  The run folders in the folder named by R012_RUNS:
#          "2021-02-07 ... (RERUN 0216)"/TABLE<j>{FDR,TPR}/N <n> {FDR,TPR} Facet <j> .xlsx
#          "2021-02-19 ... (EXCEL 0220)"/EXCEL FILES 0219/{Final Summary,FDR Facet} <j> .xlsx
#          "2021-02-24 ... (RERUN 0225)"/Final Summary <j> .xlsx
#          "2021-02-25 ... (RERUN 0309)"/EXCEL 0227/H Final Summary <j> .xlsx, EXCEL 0228/FDR Facet <j> .xlsx
#          "2021-02-17 Figure S1 plots (RERUN 0217)"/{FDR,TPR}/Table<j>{FDR,TPR}.xlsx
#          "2021-05-09 Figure S3 table"/Heatmaps AOAS.xlsx (sheet "Full Format")
#          and the published figures in the folder named by R012_FIGURES.
# Outputs: archive/simulation_summaries.csv, archive/figS1_plot_input.csv,
#          archive/figS3_full_format.csv, archive/figures/Fig*.pdf, FigS3*.png
# Usage:   R012_RUNS=/path/to/run/folders [R012_FIGURES=/path/to/figures] \
#            Rscript archive/reduce_archive.R
#          Needs the R package openxlsx.
# -----------------------------------------------------------------------------

source("R/setup.R")
library(openxlsx)

runs_dir <- Sys.getenv("R012_RUNS")
if (!nzchar(runs_dir))
  stop("Set R012_RUNS to the folder of the original runs (not part of this repository).")
figures_dir <- Sys.getenv("R012_FIGURES", file.path(dirname(dirname(runs_dir)), "5 Figures"))
RUN_FOLDERS <- c(run0207 = "2021-02-07 Figures S1 and S3 simulations (RERUN 0216)",
                 run0219 = "2021-02-19 Table 1 (C = 5) and Figure S3 simulations (EXCEL 0220)",
                 run0224 = "2021-02-24 Table 1 (C = 10) simulations (RERUN 0225)",
                 run0225 = "2021-02-25 Table 2 and Figure S3 simulations (RERUN 0309)",
                 plots0217 = "2021-02-17 Figure S1 plots (RERUN 0217)",
                 table0509 = "2021-05-09 Figure S3 table")

#' Cells of a sheet as text; numbers keep the digits a 15 digit print shows.
read_sheet <- function(..., sheet = 1) {
  cells <- read.xlsx(file.path(runs_dir, ...), sheet = sheet, colNames = FALSE,
                     skipEmptyRows = FALSE, skipEmptyCols = FALSE)
  as.data.frame(lapply(cells, function(column) {
    number <- suppressWarnings(as.numeric(column))
    ifelse(is.na(number), column, format(number, digits = 15, trim = TRUE, scientific = FALSE,
                                         drop0trailing = TRUE))
  }))
}
plain <- function(text) sub("^\\((.*)\\)$", "\\1", text)  # "(3.01)" -> "3.01"
COLUMNS <- c("study", "run", "setting", "N", "method", "gamma", "R", "R_sd", "FDR", "FDR_sd",
             "TPR", "TPR_sd")
row <- function(...) {
  given <- list(...)
  as.data.frame(lapply(setNames(COLUMNS, COLUMNS), function(name) {
    if (is.null(given[[name]]) || is.na(given[[name]])) "" else as.character(given[[name]])
  }))
}
n_of <- function(run, setting) {
  all <- simulation_settings()
  all$N[all$run == run & all$setting == setting]
}

# "Final Summary j": mean and (standard deviation) of R, FDR, TPR by method.
final_summary <- function(run, setting, ...) {
  cells <- read_sheet(...)
  stopifnot(cells[c(2, 4, 6, 8, 10, 12), 1] == c("NP", "Par", "SP", "EB TS", "EB OS", "Storey"))
  do.call(rbind, lapply(1:6, function(m) {
    mean <- cells[2 * m, 2:4]
    sd <- plain(unlist(cells[2 * m + 1, 2:4]))
    row(study = "methods", run = run, setting = setting, N = n_of(run, setting),
        method = METHODS[m], gamma = DEFAULT_GAMMA[METHODS[m]], R = mean[[1]], R_sd = sd[1],
        FDR = mean[[2]], FDR_sd = sd[2], TPR = mean[[3]], TPR_sd = sd[3])
  }))
}
# "FDR Facet j" / "TPR Facet j": one row of six means per N.
facet <- function(run, setting, N, measure, ...) {
  cells <- read_sheet(...)
  stopifnot(cells[1, ] == c("Non-parametric", "Parametric", "Semi-parametric", "Two Stage",
                            "One Stage", "Storey's"), nrow(cells) == length(N) + 1)
  do.call(rbind, lapply(seq_along(N), function(k) {
    out <- do.call(rbind, lapply(1:6, function(m) {
      row(study = "methods", run = run, setting = setting, N = N[k], method = METHODS[m],
          gamma = DEFAULT_GAMMA[METHODS[m]])
    }))
    out[[measure]] <- unlist(cells[k + 1, ])
    out
  }))
}
merge_measures <- function(a, b) {  # same rows, different measures filled
  stopifnot(a[1:6] == b[1:6])
  for (name in COLUMNS[7:12]) a[[name]] <- ifelse(a[[name]] == "", b[[name]], a[[name]])
  a
}

# run0219 (C = 5): the FDR of settings 1 to 30, and the full summary of the
# five settings of Table 1.
run0219 <- do.call(rbind, lapply(1:30, function(j) {
  fdr <- facet("run0219", j, n_of("run0219", j), "FDR", RUN_FOLDERS[["run0219"]],
               "EXCEL FILES 0219", sprintf("FDR Facet %d .xlsx", j))
  if (!j %in% c(23, 5, 29, 11, 17)) return(fdr)
  full <- final_summary("run0219", j, RUN_FOLDERS[["run0219"]], "EXCEL FILES 0219",
                        sprintf("Final Summary %d .xlsx", j))
  stopifnot(as.numeric(full$FDR) == as.numeric(fdr$FDR))
  full
}))
# run0224 (C = 10): the five settings of Table 1.
run0224 <- do.call(rbind, lapply(c(1, 2, 3, 6, 11), function(j) {
  final_summary("run0224", j, RUN_FOLDERS[["run0224"]], sprintf("Final Summary %d .xlsx", j))
}))
# Prior sensitivity (Table 2): settings of run0224, NP and SP for five priors.
prior <- do.call(rbind, lapply(c(1, 2, 3, 6, 11), function(j) {
  cells <- read_sheet(RUN_FOLDERS[["run0225"]], "EXCEL 0227", sprintf("H Final Summary %d .xlsx", j))
  do.call(rbind, lapply(seq_along(PRIOR_GAMMAS), function(k) {
    do.call(rbind, lapply(1:2, function(m) {
      mean <- unlist(cells[2 * k, 3 * m - 1 + 0:2])
      sd <- plain(unlist(cells[2 * k + 1, 3 * m - 1 + 0:2]))
      row(study = "prior", run = "run0224", setting = j, N = n_of("run0224", j),
          method = c("NP", "SP")[m], gamma = PRIOR_GAMMAS[[k]], R = mean[1], R_sd = sd[1],
          FDR = mean[2], FDR_sd = sd[2], TPR = mean[3], TPR_sd = sd[3])
    }))
  }))
}))
# run0227 (C = 5): the FDR of settings 1 to 20.
run0227 <- do.call(rbind, lapply(1:20, function(j) {
  facet("run0227", j, n_of("run0227", j), "FDR", RUN_FOLDERS[["run0225"]], "EXCEL 0228",
        sprintf("FDR Facet %d .xlsx", j))
}))
# run0207: FDR and TPR of settings 16 and 18, FDR of settings 10 and 14. Each
# file holds the N of one call of the run.
blocks <- list("50" = c(50, 100, 200), "300" = c(300, 400), "500" = c(500, 1000),
               "600" = c(600, 700), "800" = c(800, 900))
run0207_measure <- function(setting, measure) {
  out <- do.call(rbind, lapply(names(blocks), function(first) {
    facet("run0207", setting, blocks[[first]], measure, RUN_FOLDERS[["run0207"]],
          sprintf("TABLE%d%s", setting, measure), sprintf("N %s %s Facet %d .xlsx", first, measure, setting))
  }))
  out <- out[as.numeric(out$N) %in% n_of("run0207", setting), ]
  out[order(as.numeric(out$N)), ]
}
run0207 <- rbind(run0207_measure(10, "FDR"), run0207_measure(14, "FDR"),
                 merge_measures(run0207_measure(16, "FDR"), run0207_measure(16, "TPR")),
                 merge_measures(run0207_measure(18, "FDR"), run0207_measure(18, "TPR")))

write.csv(rbind(run0207, run0219, run0224, prior, run0227),
          file.path(PATHS$archive, "simulation_summaries.csv"), row.names = FALSE, quote = FALSE)

# Figure S1: the four tables the plots were drawn from.
plot_input <- do.call(rbind, lapply(c(16, 18), function(setting) {
  do.call(rbind, lapply(c("FDR", "TPR"), function(measure) {
    cells <- read_sheet(RUN_FOLDERS[["plots0217"]], measure, sprintf("Table%d%s.xlsx", setting, measure))
    data.frame(setting = setting, measure = measure,
               setNames(cells[1:11, ], c("C", "N", METHODS)))
  }))
}))
write.csv(plot_input, file.path(PATHS$archive, "figS1_plot_input.csv"), row.names = FALSE, quote = FALSE)

# Figure S3: the sheet the published table shows. Rows N x lambda; for each
# theta the columns SP, P, NP, OS, Storey, TS.
sheet <- read_sheet(RUN_FOLDERS[["table0509"]], "Heatmaps AOAS.xlsx", sheet = "Full Format")
stopifnot(sheet[2, 3:8] == c("SP", "P", "NP", "OS", "St", "TS"))
full_format <- do.call(rbind, lapply(3:22, function(r) {
  do.call(rbind, lapply(1:3, function(k) {
    data.frame(N = sheet[r - (r - 3) %% 4, 1], lambda = sheet[r, 2], theta = c(0.15, 0.2, 0.3)[k],
               method = c("SP", "GP", "NP", "OS", "Storey", "TS"),
               FDR = unlist(sheet[r, 2 + 6 * (k - 1) + 1:6]))
  }))
}))
write.csv(full_format, file.path(PATHS$archive, "figS3_full_format.csv"), row.names = FALSE, quote = FALSE)

# Published figure files, under the names the scripts give their figures.
published <- list.files(figures_dir, recursive = TRUE)
copy_figure <- function(pattern, name) {
  from <- published[grepl(pattern, published, fixed = TRUE)]
  stopifnot(length(from) == 1)
  dir.create(file.path(PATHS$archive, "figures"), showWarnings = FALSE)
  file.copy(file.path(figures_dir, from), file.path(PATHS$archive, "figures", name), overwrite = TRUE,
            copy.mode = FALSE)
}
for (i in seq_along(DOMAINS)) {
  copy_figure(sprintf("1%s Bar plot", LETTERS[i]), fig1_file(i))
}
for (i in 1:4) {
  copy_figure(sprintf("9 Supplement/1%s ", LETTERS[i]), figS1_file(i))
}
for (i in 1:9) {
  copy_figure(sprintf("9 Supplement/2%s ", LETTERS[i]), figS2_file(i))
}
copy_figure("9 Supplement/3 FDR", "FigS3_fdr_table.png")
