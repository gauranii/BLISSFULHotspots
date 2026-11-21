# -----------------------------------------------------------------------------
# File:    tests/testthat/test-simulate.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Unit tests of R/simulate.R: the simulated mixture under a seed and
#          against its parameters, the settings of the runs, one replication
#          of each method, the stored replications and their summaries on
#          files made by hand, a small study end to end, the archived
#          summaries and the layout of Tables 1 and 2.
# Inputs:  archive/simulation_summaries.csv
# Outputs: None (replications are written to temporary folders).
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/setup.R")

#' Write the replication files of a cell from a data frame, one row each.
write_replications <- function(cell, replications, results) {
  folder <- replication_dir(cell, cell$method, cell$gamma, SMALL_SETTINGS, results)
  dir.create(folder, recursive = TRUE)
  for (d in seq_len(nrow(replications))) {
    write.csv(replications[d, ], file.path(folder, sprintf("rep%04d.csv", d)), row.names = FALSE)
  }
  folder
}

test_that("rmixture gives the same counts under the same seed, whatever came before", {
  first <- rmixture(80, 0.4, 2, 0.3, 0.8, "binomial", C = 5, seed = 1440)
  runif(3)
  second <- rmixture(80, 0.4, 2, 0.3, 0.8, "binomial", C = 5, seed = 1440)
  expect_identical(first, second)
  expect_false(identical(first$draw, rmixture(80, 0.4, 2, 0.3, 0.8, "binomial", C = 5, seed = 1441)$draw))
  expect_named(first, c("draw", "true_grp"))
  expect_length(first$draw, 80)
  expect_true(all(first$true_grp %in% c(0, 1)))
})

test_that("rmixture has the null proportion and the null mean it is given", {
  mixture <- rmixture(5000, 0.4, 2, 0.3, 0.8, "binomial", C = 5, seed = 3)
  expect_equal(mean(mixture$true_grp == 0), 0.8, tolerance = 0.03)
  null <- mixture$draw[mixture$true_grp == 0]
  expect_equal(mean(null), 0.6 * 2 / 0.7, tolerance = 0.05)
  expect_equal(mean(null == 0), 0.4 + 0.6 * exp(-2), tolerance = 0.05)
  # the non-null counts lie above the cutoff: a ZIGP count plus Binomial(250, 0.2) + C - 39
  alt <- mixture$draw[mixture$true_grp == 1]
  expect_true(all(alt > 5))
  expect_equal(mean(alt), 0.6 * 2 / 0.7 + 50 + 5 - 39, tolerance = 0.03)
  # all null, all non-null
  expect_equal(rmixture(30, 0.4, 2, 0.3, 1, "binomial", 5, seed = 1)$true_grp, rep(0, 30))
  expect_equal(rmixture(30, 0.4, 2, 0.3, 0, "binomial", 5, seed = 1)$true_grp, rep(1, 30))
})

test_that("rmixture draws every f1 above the cutoff", {
  for (f1 in c("binomial", "geometric", "binomial mixture", "poisson mixture")) {
    mixture <- rmixture(400, 0.4, 2, 0.3, 0.5, f1, C = 7, seed = 9)
    expect_true(all(mixture$draw[mixture$true_grp == 1] > 7), info = f1)
    expect_equal(mean(mixture$true_grp), 0.5, tolerance = 0.15, info = f1)
  }
  # geometric f1: ZIGP + Geometric(p) + C + 1, with mean (1 - p) / p
  geometric <- rmixture(3000, 0.4, 2, 0.3, 0, "geometric", C = 7, seed = 9, geometric_prob = 0.5)
  expect_equal(mean(geometric$draw), 0.6 * 2 / 0.7 + 1 + 8, tolerance = 0.03)
  # Poisson mixture f1: mean 0.4 * 5 + 0.6 * 10
  poisson <- rmixture(3000, 0.4, 2, 0.3, 0, "poisson mixture", C = 7, seed = 9)
  expect_equal(mean(poisson$draw), 0.6 * 2 / 0.7 + 8 + 8, tolerance = 0.03)
  expect_error(rmixture(10, 0.4, 2, 0.3, 0.5, "uniform", 5, seed = 1), "unknown f1")
  expect_error(rmixture(10, 1.4, 2, 0.3, 0.5, "binomial", 5, seed = 1))
})

test_that("simulated_data splits the counts into the two groups", {
  setting <- small_setting(N = 80, lambda = 2, theta = 0.3, C = 5)
  data <- simulated_data(setting, seed = 1440)
  mixture <- rmixture(80, 0.4, 2, 0.3, 0.8, "binomial", C = 5, seed = 1440)
  expect_named(data, c("null", "alt", "x", "values", "null_frequency", "alt_frequency"))
  expect_equal(data$null, mixture$draw[mixture$true_grp == 0])
  expect_equal(data$alt, mixture$draw[mixture$true_grp == 1])
  expect_equal(data$x, c(data$null, data$alt))
  expect_equal(data$values, sort(unique(mixture$draw)))
  expect_equal(sum(data$null_frequency), length(data$null))
  expect_equal(sum(data$alt_frequency), length(data$alt))
  expect_equal(data$null_frequency + data$alt_frequency, as.numeric(table(mixture$draw)))
  expect_identical(simulated_data(setting, seed = 1440), data)
  expect_error(simulated_data(transform(setting, pi0 = 1), seed = 1), "fewer than two counts")
})

test_that("simulation_settings lists the settings of the four runs", {
  all <- simulation_settings()
  expect_named(all, c("run", "setting", "eta", "lambda", "theta", "pi0", "f1", "C", "N"))
  expect_equal(as.vector(table(all$run)[names(SIMULATION_RUNS)]), c(32, 30, 11, 20))
  expect_true(all(all$eta == 0.4 & all$f1 == "binomial"))
  expect_equal(sort(all$setting[all$run == "run0219"]), 1:30)
  expect_equal(sort(all$setting[all$run == "run0227"]), 1:20)
  expect_equal(all$N[all$run == "run0224"], c(50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000))
  expect_equal(unique(all$pi0[all$setting == 16 & all$run == "run0207"]), 0.35)
  expect_equal(unique(all$pi0[!(all$setting == 16 & all$run == "run0207")]), 0.8)
  # no setting is listed twice
  expect_equal(anyDuplicated(all[c("run", "setting", "N")]), 0)
  expect_equal(unlist(all[all$run == "run0219" & all$setting == 5, c("lambda", "theta", "C", "N")],
                      use.names = FALSE), c(2, 0.3, 5, 100))
})

test_that("SIMULATION_RUNS and the method constants name what the settings use", {
  expect_true(all(vapply(SIMULATION_RUNS, function(run) run$version, "") %in% names(SAMPLER_VERSIONS)))
  expect_equal(SIMULATION_RUNS$run0224$bayes_seed, 1440)
  expect_equal(METHODS, c("NP", "GP", "SP", "TS", "OS", "Storey"))
  expect_named(METHOD_LABELS, METHODS)
  expect_equal(unname(PRIOR_GAMMAS), c(0.02, 0.5, 1, 1.5, 3))
  expect_equal(DEFAULT_GAMMA, c(NP = 1, GP = 0.5, SP = 1))
})

test_that("item_settings picks the settings of each table and figure", {
  table1 <- item_settings("table1")
  expect_equal(nrow(table1), 10)
  expect_equal(sort(table1$N[table1$C == 10]), c(50, 100, 200, 500, 1000))
  expect_equal(sort(table1$N[table1$C == 5]), c(50, 100, 200, 500, 1000))
  expect_true(all(table1$lambda == 2 & table1$theta == 0.3))
  table2 <- item_settings("table2")
  expect_equal(nrow(table2), 5)
  expect_true(all(table2$run == "run0224"))
  figS1 <- item_settings("figS1")
  expect_equal(nrow(figS1), 22)
  expect_equal(sort(unique(figS1$setting)), c(16, 18))
  expect_equal(nrow(item_settings("figS3")), 60)
  expect_error(item_settings("table9"), "unknown item")
})

test_that("bayes_replication counts the rejections in the data set of its seed", {
  setting <- small_setting(N = 60)
  first <- bayes_replication(setting, "NP", 1, d = 3, SMALL_SETTINGS)
  runif(2)
  expect_identical(bayes_replication(setting, "NP", 1, d = 3, SMALL_SETTINGS), first)
  expect_named(first, c("d", "R", "V", "TP", "n_alt"))
  data <- simulated_data(setting, SIMULATION_RUNS$run0224$bayes_seed + 3)
  expect_equal(first$d, 3)
  expect_equal(first$n_alt, length(data$alt))
  expect_equal(first$R, first$V + first$TP)
  expect_equal(c(first$V, first$TP), c(0, length(data$alt)))   # the groups are far apart
  # the chain continues the random numbers of the data set
  fit <- bayes_local_fdr(simulated_data(setting, 1440 + 3)$x, "NP", 1, sampler_settings(2, 0.25),
                         draws = 30, burn = 10, version = "2021-02-19")
  expect_equal(first$R, sum(fit$frequency[fit$fdr <= 0.05]))
  expect_equal(bayes_replication(setting, "NP", 1, d = 4, SMALL_SETTINGS)$d, 4)
})

test_that("eb_replication counts the rejections of the three procedures", {
  setting <- small_setting(N = 400)
  first <- eb_replication(setting, d = 2, SMALL_SETTINGS)
  expect_named(first, c("d", "n_alt", "TS_R", "TS_V", "TS_TP", "OS_R", "OS_V", "OS_TP",
                        "Storey_R", "Storey_V", "Storey_TP"))
  data <- simulated_data(setting, SIMULATION_RUNS$run0224$eb_seed + 2)
  expect_equal(first$d, 2)
  expect_equal(first$n_alt, length(data$alt))
  expect_equal(first[-(1:2)], eb_simulation_counts(data, 10, 0.05))
  expect_identical(eb_replication(setting, d = 2, SMALL_SETTINGS), first)
})

test_that("replication_dir names the data, the method and what the chain depends on", {
  setting <- small_setting(N = 400)
  expect_equal(replication_dir(setting, "NP", 1, PAPER, "res"),
               paste0("res/replications/eta0.4_lambda1.5_theta0_pi00.8_binomial_C10_N400/",
                      "NP_gamma1_seed1440_sampler2021-02-19_init0.25_draws5000_burn1000"))
  expect_equal(replication_dir(setting, "EB", NA, PAPER, "res"),
               "res/replications/eta0.4_lambda1.5_theta0_pi00.8_binomial_C10_N400/EB_seed889191")
  # shorter chains and another run are stored elsewhere
  expect_match(replication_dir(setting, "SP", 0.02, SMALL_SETTINGS, "res"),
               "SP_gamma0.02_seed1440_sampler2021-02-19_init0.25_draws30_burn10$")
  expect_match(replication_dir(transform(setting, run = "run0207"), "GP", 0.5, PAPER, "res"),
               "GP_gamma0.5_seed889191_sampler2021-02-07_init0.2_draws5000_burn1000$")
})

test_that("run_replications writes one file per replication and keeps the finished ones", {
  results <- new_folder()
  cells <- rbind(cbind(small_setting(N = 60), method = "NP", gamma = 1),
                 cbind(small_setting(N = 60), method = "EB", gamma = NA))
  expect_output(run_replications(cells, 1:2, SMALL_SETTINGS, results), "4 replications to run, 1 cores")
  bayes <- file.path(replication_dir(cells[1, ], "NP", 1, SMALL_SETTINGS, results), "rep0002.csv")
  eb <- file.path(replication_dir(cells[2, ], "EB", NA, SMALL_SETTINGS, results), "rep0001.csv")
  expect_true(all(file.exists(bayes, eb)))
  expect_equal(read.csv(bayes), bayes_replication(cells[1, ], "NP", 1, 2, SMALL_SETTINGS))
  expect_equal(read.csv(eb), eb_replication(cells[2, ], 1, SMALL_SETTINGS))
  # a second call runs only what is missing, and gives the same files in parallel processes
  writeLines("kept", bayes)
  expect_output(run_replications(cells, 1:3, SMALL_SETTINGS, results, cores = 2),
                "2 replications to run, 2 cores")
  expect_equal(readLines(bayes), "kept")
  third <- file.path(dirname(bayes), "rep0003.csv")
  expect_equal(read.csv(third), bayes_replication(cells[1, ], "NP", 1, 3, SMALL_SETTINGS))
})

test_that("run_replications stops when a replication fails", {
  cells <- cbind(small_setting(N = 60, pi0 = 1), method = "EB", gamma = NA)   # no non-null count
  expect_error(suppressWarnings(quietly(run_replications(cells, 1:2, SMALL_SETTINGS, new_folder(),
                                                        cores = 2))),
               "2 replications failed.*fewer than two counts")
})

test_that("summarize_replications averages a Bayesian cell worked by hand", {
  results <- new_folder()
  cell <- cbind(small_setting(), method = "NP", gamma = 1)
  write_replications(cell, data.frame(d = 1:3, R = c(4, 0, 10), V = c(1, 0, 0), TP = c(3, 0, 10),
                                      n_alt = c(6, 5, 10)), results)
  summary <- summarize_replications(cell, 3, SMALL_SETTINGS, results)
  # R 4, 0, 10; FDR 1/4, 0 (no rejection), 0; TPR 3/6, 0, 1
  expect_equal(summary, data.frame(run = "run0224", setting = 1, N = 60, method = "NP", gamma = 1,
                                   R = 4.67, R_sd = 5.03, FDR = 0.0833, FDR_sd = 0.1443,
                                   TPR = 0.5, TPR_sd = 0.5))
  # the first two replications only
  expect_equal(summarize_replications(cell, 2, SMALL_SETTINGS, results)$R, 2)
  expect_error(summarize_replications(cell, 4, SMALL_SETTINGS, results), "replications missing")
})

test_that("summarize_replications gives one row per procedure of an EB cell", {
  results <- new_folder()
  cell <- cbind(small_setting(), method = "EB", gamma = NA)
  write_replications(cell, data.frame(d = 1:2, n_alt = c(10, 10),
                                      TS_R = c(10, 0), TS_V = c(2, 0), TS_TP = c(8, 0),
                                      OS_R = c(5, 5), OS_V = c(0, 1), OS_TP = c(5, 4),
                                      Storey_R = c(8, 0), Storey_V = c(2, 0), Storey_TP = c(6, 0)), results)
  summary <- summarize_replications(cell, 2, SMALL_SETTINGS, results)
  expect_equal(summary$method, c("TS", "OS", "Storey"))
  expect_true(all(is.na(summary$gamma)))
  expect_equal(summary$R, c(5, 5, 4))
  expect_equal(summary$FDR, c(0.1, 0.1, 0.25))   # Storey: the replication without rejections is left out
  expect_equal(summary$TPR, c(0.4, 0.45, 0.3))
  expect_equal(summary$FDR_sd, c(round(sd(c(0.2, 0)), 4), round(sd(c(0, 0.2)), 4), NA))
})

test_that("simulation_summaries runs a small study and summarizes it when complete", {
  paths <- list(results = new_folder())
  opts <- list(rerun = TRUE, runs = "1", replications = 2, cores = 2)
  expect_output(first <- simulation_summaries("table1", opts, SMALL_SETTINGS, paths),
                "replications 1 to 2 are not all stored yet")
  expect_null(first)
  summaries <- quietly(simulation_summaries("table1", modifyList(opts, list(runs = "2")),
                                            SMALL_SETTINGS, paths))
  # 10 settings: three Bayesian methods and the three procedures of the EB fit
  expect_equal(nrow(summaries), 60)
  expect_named(summaries, c("run", "setting", "N", "method", "gamma", "R", "R_sd", "FDR", "FDR_sd",
                            "TPR", "TPR_sd"))
  expect_equal(as.vector(table(summaries$method)[METHODS]), rep(10, 6))
  expect_true(all(summaries$FDR >= 0 & summaries$FDR <= 1, na.rm = TRUE))
  expect_true(all(summaries$TPR >= 0 & summaries$TPR <= 1))
  expect_equal(dim(table1_layout(summaries)), c(30L, 14L))
  # nothing is left to run, and the summaries are read again from the files
  expect_output(again <- simulation_summaries("table1", modifyList(opts, list(runs = "2")),
                                              SMALL_SETTINGS, paths), "0 replications to run")
  expect_identical(again, summaries)
  one <- summaries[summaries$run == "run0224" & summaries$N == 50 & summaries$method == "SP", ]
  setup <- item_settings("table1")
  cell <- cbind(setup[setup$run == "run0224" & setup$N == 50, ], method = "SP", gamma = 1)
  rownames(one) <- NULL
  expect_equal(one, summarize_replications(cell, 2, SMALL_SETTINGS, paths$results))
})

test_that("archived_summaries reads the runs behind each table and figure", {
  table1 <- archived_summaries("table1", PATHS)
  expect_equal(nrow(table1), 60)
  expect_equal(as.vector(table(table1$method)[METHODS]), rep(10, 6))
  expect_true(all(table1$study == "methods"))
  table2 <- archived_summaries("table2", PATHS)
  expect_equal(nrow(table2), 50)
  expect_true(all(table2$study == "prior" & table2$method %in% c("NP", "SP")))
  expect_equal(sort(unique(table2$gamma)), unname(PRIOR_GAMMAS))
  expect_equal(nrow(archived_summaries("figS1", PATHS)), 132)   # 22 settings, 6 methods
  expect_equal(nrow(archived_summaries("figS3", PATHS)), 360)   # 60 settings, 6 methods
})

test_that("with_parameters adds the parameters of each setting", {
  summaries <- data.frame(run = c("run0224", "run0219", "run0207"), setting = c(1, 5, 16),
                          N = c(50, 100, 300), FDR = c(0.01, 0.02, 0.03))
  out <- with_parameters(summaries)
  expect_named(out, c("run", "setting", "N", "FDR", "lambda", "theta", "pi0", "C"))
  expect_equal(out$lambda, c(2, 2, 1.5))
  expect_equal(out$theta, c(0.3, 0.3, 0))
  expect_equal(out$pi0, c(0.8, 0.8, 0.35))
  expect_equal(out$C, c(10, 5, 10))
})

test_that("summary_block takes the six summary columns and prefixes them", {
  summaries <- data.frame(N = 50, R = 9.33, R_sd = 3.01, FDR = 0.0141, FDR_sd = 0.0496,
                          TPR = 0.9404, TPR_sd = 0.0919, method = "NP")
  expect_equal(summary_block(summaries, "C10"),
               data.frame(C10_R = 9.33, C10_R_sd = 3.01, C10_FDR = 0.0141, C10_FDR_sd = 0.0496,
                          C10_TPR = 0.9404, C10_TPR_sd = 0.0919))
})

test_that("table1_layout has one row per N and method, with C = 10 beside C = 5", {
  archive <- archived_summaries("table1", PATHS)
  table <- table1_layout(archive)
  expect_equal(dim(table), c(30L, 14L))
  expect_equal(table$N, rep(c(50, 100, 200, 500, 1000), each = 6))
  expect_equal(table$method, rep(unname(METHOD_LABELS), 5))
  # a cell of each block against the archive
  c10 <- archive[archive$run == "run0224" & archive$N == 200 & archive$method == "TS", ]
  c5 <- archive[archive$run == "run0219" & archive$N == 200 & archive$method == "TS", ]
  row <- table[table$N == 200 & table$method == "Two-Stage", ]
  expect_equal(c(row$C10_R, row$C10_FDR, row$C10_TPR_sd), c(c10$R, c10$FDR, c10$TPR_sd))
  expect_equal(c(row$C5_R, row$C5_FDR, row$C5_TPR_sd), c(c5$R, c5$FDR, c5$TPR_sd))
  expect_equal(unlist(table[1, c("C10_R", "C10_FDR", "C10_TPR")], use.names = FALSE),
               c(9.33, 0.0141, 0.9404))
})

test_that("table2_layout has one row per N and prior, with NP beside SP", {
  archive <- archived_summaries("table2", PATHS)
  table <- table2_layout(archive)
  expect_equal(dim(table), c(25L, 15L))
  expect_equal(table$N, rep(c(50, 100, 200, 500, 1000), each = 5))
  expect_equal(table$prior, rep(names(PRIOR_GAMMAS), 5))
  expect_equal(table$gamma, rep(unname(PRIOR_GAMMAS), 5))
  np <- archive[archive$N == 500 & archive$gamma == 3 & archive$method == "NP", ]
  sp <- archive[archive$N == 500 & archive$gamma == 3 & archive$method == "SP", ]
  row <- table[table$N == 500 & table$prior == "Informative", ]
  expect_equal(c(row$NP_R, row$NP_FDR, row$SP_R, row$SP_TPR), c(np$R, np$FDR, sp$R, sp$TPR))
  expect_equal(unlist(table[1, c("NP_R", "NP_FDR", "SP_TPR")], use.names = FALSE), c(10.70, 0.0718, 0.7699))
})
