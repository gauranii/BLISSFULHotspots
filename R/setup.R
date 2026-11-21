# -----------------------------------------------------------------------------
# File:    R/setup.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. All figures and tables.
# Purpose: Loader sourced first by every script: finds the repository root,
#          sources the library files in R/, and defines the settings of the
#          paper, the project paths and the command line parser.
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R")   # from a script started in the repository root
# Source:  sim_initial0225.R, fcns6_initial.R, real_scenarios.R (settings of
#          the original runs).
# -----------------------------------------------------------------------------

#' File being read by `source()`.
#'
#' @return The path given to the innermost running `source()` call.
sourced_file <- function() {
  for (frame in rev(sys.frames())) {
    if (!is.null(frame$ofile)) return(frame$ofile)
  }
  stop("R/setup.R must be loaded with source()")
}

REPO_ROOT <- normalizePath(file.path(dirname(sourced_file()), ".."))

for (file in c("distributions.R", "bayes_fdr.R", "eb_fdr.R", "simulate.R", "real_data.R",
               "plots.R")) {
  source(file.path(REPO_ROOT, "R", file))
}

# ----- Settings of the paper, passed to the functions as `settings` -----
PAPER <- list(
  alpha        = 0.05,   # level of the local fdr and of the FDR
  replications = 1000,   # simulated data sets per setting
  sim_draws    = 5000,   # stored iterations of the sampler, simulations
  sim_burn     = 1000,   # burn-in iterations, simulations
  real_draws   = 10000,  # stored iterations, protein domains
  real_burn    = 1000,   # burn-in iterations, protein domains
  real_seed    = 1234,   # seed set before each model, protein domains
  # Dirichlet concentration (NP, SP) per model in the protein domain analysis
  real_gamma   = c(NP = 1.5, SP = 0.02, P = 0.5, GP = 0.5),
  figS2_seed   = 1439,   # seed of the simulated data sets of Figure S2
  figure_size  = 7       # width and height of every figure, inches
)

# ----- Paths -----

#' Folders and files of the project.
#'
#' @param root Repository root.
#' @param output Folder in which results/ and figures/ are written.
#' @return Named list of paths: `domains`, `archive`, `results`, `figures`.
project_paths <- function(root, output = Sys.getenv("LOCALFDR_OUTPUT", root)) {
  list(domains = file.path(root, "data", "domains"),
       archive = file.path(root, "archive"),
       results = file.path(output, "results"),
       figures = file.path(output, "figures"))
}

PATHS <- project_paths(REPO_ROOT)

#' Path of an output file, creating its folder when needed.
#'
#' @param folder Output folder (for example `paths$results`).
#' @param name File name.
#' @return The path `folder/name`.
output_file <- function(folder, name) {
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  file.path(folder, name)
}

# ----- Command line -----

#' Read `--name=value` arguments and bare `--name` flags.
#'
#' Every script accepts `--rerun`: recompute by simulation instead of reading
#' the archived runs.
#'
#' @param defaults Named list of options and their default values. A value
#'   given on the command line is converted to the type of its default; a
#'   flag sets a logical default to TRUE.
#' @param args Command line arguments.
#' @return `defaults` with the values given on the command line, and `rerun`.
parse_args <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  defaults <- c(list(rerun = FALSE), defaults)
  pairs <- regmatches(args, regexec("^--([a-z_]+)(=(.*))?$", args))
  if (any(lengths(pairs) != 4)) stop("arguments must look like --name=value or --name")
  for (pair in pairs) {
    name <- pair[2]
    if (!name %in% names(defaults)) stop("unknown argument --", name)
    is_flag <- pair[3] == ""
    if (is_flag != is.logical(defaults[[name]])) {
      stop("--", name, if (is_flag) " needs a value" else " takes no value")
    }
    defaults[[name]] <- if (is_flag) TRUE else methods::as(pair[4], class(defaults[[name]]))
  }
  defaults
}

#' Default options of a script that can rerun a simulation study.
#'
#' @param settings Settings of the paper (`replications`).
#' @return Named list of defaults for `parse_args()`: `runs` (replications to
#'   run now), `replications` (replications of the study) and `cores`.
simulation_options <- function(settings) {
  list(runs = as.character(settings$replications), replications = settings$replications, cores = 4)
}

#' Turn an index text into indices.
#'
#' @param text "3:7" (3 to 7), "2,5,9", or a single number N.
#' @param count If TRUE a single number N means 1 to N.
#' @return Integer vector.
parse_index_set <- function(text, count = FALSE) {
  if (grepl("^[0-9]+:[0-9]+$", text)) {
    ends <- as.integer(strsplit(text, ":")[[1]])
    return(ends[1]:ends[2])
  }
  index <- as.integer(strsplit(text, ",")[[1]])
  if (count && length(index) == 1) seq_len(index) else index
}

#' Name of an output of a simulation study.
#'
#' The name carries the number of replications unless it is that of the paper.
#'
#' @param name Name of the table or figure output, without extension.
#' @param opts Options with `rerun` and `replications`.
#' @param settings Settings of the paper.
#' @return `name`, or `name_first<D>-replications`.
study_name <- function(name, opts, settings) {
  if (!opts$rerun || opts$replications == settings$replications) return(name)
  sprintf("%s_first%d-replications", name, opts$replications)
}

#' Print a time-stamped message.
#'
#' @param msg Format of the message, as for `sprintf()`.
#' @param ... Values for the format.
#' @return NULL, invisibly.
log_message <- function(msg, ...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", sprintf(msg, ...), "\n")
}
