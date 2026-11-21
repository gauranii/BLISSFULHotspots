# -----------------------------------------------------------------------------
# File:    tests/testthat/test-setup.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Unit tests of R/setup.R: the loader, the settings, the project
#          paths, the command line parser, the index and name helpers and
#          the log line.
# Inputs:  None (files are written to temporary folders).
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/setup.R")

test_that("sourced_file is the file of the innermost source() call", {
  folder <- new_folder()
  writeLines("inner <- basename(sourced_file())", file.path(folder, "inner.R"))
  writeLines(c(sprintf("source('%s', local = TRUE)", file.path(folder, "inner.R")),
               "outer <- basename(sourced_file())"), file.path(folder, "outer.R"))
  seen <- new.env()
  source(file.path(folder, "outer.R"), local = seen)
  expect_equal(c(seen$inner, seen$outer), c("inner.R", "outer.R"))
  expect_error(sourced_file(), "must be loaded with source")
})

test_that("setup.R finds the repository and loads every module", {
  expect_true(file.exists(file.path(REPO_ROOT, "flake.nix")))
  for (name in c("dzigp", "bayes_local_fdr", "eb_fit", "rmixture", "read_domain", "rate_plot")) {
    expect_true(is.function(get(name)), info = name)
  }
})

test_that("PAPER holds the settings of the paper", {
  expect_equal(PAPER$alpha, 0.05)
  expect_equal(PAPER$replications, 1000)
  expect_equal(c(PAPER$sim_draws, PAPER$sim_burn), c(5000, 1000))
  expect_equal(c(PAPER$real_draws, PAPER$real_burn, PAPER$real_seed), c(10000, 1000, 1234))
  expect_equal(PAPER$real_gamma, c(NP = 1.5, SP = 0.02, P = 0.5, GP = 0.5))
  expect_equal(PAPER$figS2_seed, 1439)
})

test_that("project_paths puts the inputs under the root and the outputs under output", {
  paths <- project_paths("/repo", "/out")
  expect_equal(paths, list(domains = "/repo/data/domains", archive = "/repo/archive",
                           results = "/out/results", figures = "/out/figures"))
  expect_equal(names(PATHS), names(paths))
  expect_true(dir.exists(PATHS$domains))
})

test_that("project_paths reads LOCALFDR_OUTPUT, and the root without it", {
  old <- Sys.getenv("LOCALFDR_OUTPUT", NA)
  on.exit(if (is.na(old)) Sys.unsetenv("LOCALFDR_OUTPUT") else Sys.setenv(LOCALFDR_OUTPUT = old))
  Sys.setenv(LOCALFDR_OUTPUT = "/elsewhere")
  expect_equal(project_paths("/repo")$results, "/elsewhere/results")
  Sys.unsetenv("LOCALFDR_OUTPUT")
  expect_equal(project_paths("/repo")$figures, "/repo/figures")
})

test_that("output_file creates the folder and returns the path", {
  folder <- file.path(new_folder(), "a", "b")
  expect_false(dir.exists(folder))
  expect_equal(output_file(folder, "x.csv"), file.path(folder, "x.csv"))
  expect_true(dir.exists(folder))
  expect_equal(output_file(folder, "y.csv"), file.path(folder, "y.csv"))  # folder exists already
})

test_that("parse_args converts values to the type of the default and sets flags", {
  defaults <- list(cores = 4, domains = "all", poisson = FALSE)
  expect_equal(parse_args(defaults, character(0)), c(list(rerun = FALSE), defaults))
  opts <- parse_args(defaults, c("--cores=2", "--domains=cd00083,cd00054", "--poisson", "--rerun"))
  expect_equal(opts, list(rerun = TRUE, cores = 2, domains = "cd00083,cd00054", poisson = TRUE))
  expect_type(opts$cores, "double")
  expect_equal(parse_args(list(runs = "1000"), "--runs=3:7")$runs, "3:7")
})

test_that("parse_args refuses what it does not know", {
  defaults <- list(cores = 4, poisson = FALSE)
  expect_error(parse_args(defaults, "--corse=2"), "unknown argument --corse")
  expect_error(parse_args(defaults, "cores=2"), "must look like")
  expect_error(parse_args(defaults, "--cores"), "--cores needs a value")
  expect_error(parse_args(defaults, "--poisson=yes"), "--poisson takes no value")
})

test_that("simulation_options are the replications of the study and 4 cores", {
  expect_equal(simulation_options(list(replications = 250)),
               list(runs = "250", replications = 250, cores = 4))
  expect_equal(parse_args(simulation_options(PAPER), c("--rerun", "--runs=1:2", "--replications=2")),
               list(rerun = TRUE, runs = "1:2", replications = 2, cores = 4))
})

test_that("parse_index_set reads a range, a list and a count", {
  expect_equal(parse_index_set("3:7"), 3:7)
  expect_equal(parse_index_set("2,5,9"), c(2L, 5L, 9L))
  expect_equal(parse_index_set("4"), 4L)
  expect_equal(parse_index_set("4", count = TRUE), 1:4)
  expect_equal(parse_index_set("4,6", count = TRUE), c(4L, 6L))
})

test_that("study_name carries the number of replications unless it is the paper's", {
  settings <- list(replications = 1000)
  expect_equal(study_name("Table1", list(rerun = FALSE, replications = 10), settings), "Table1")
  expect_equal(study_name("Table1", list(rerun = TRUE, replications = 1000), settings), "Table1")
  expect_equal(study_name("Table1", list(rerun = TRUE, replications = 10), settings),
               "Table1_first10-replications")
})

test_that("log_message prints the time and the formatted message", {
  expect_output(log_message("wrote %d files to %s", 3, "here"),
                "^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9:]{8} - wrote 3 files to here")
  expect_null(quietly(log_message("x")))
})
