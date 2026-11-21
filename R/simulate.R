# -----------------------------------------------------------------------------
# File:    R/simulate.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Tables 1 and 2, Figures S1 to S3.
# Purpose: Simulated count data from the two-group mixture, the settings of the
#          simulation runs of the paper, and the driver that runs replications
#          in parallel, stores one file per replication and summarizes them.
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R")
# Source:  simFDR_fcns.R (rmixture, bayes_fdr, EM_fdr, facet, prior) and the
#          settings files initial_vals0207.R, sim_initial0219.R,
#          sim_initial0225.R, sim_initial0227.R of the simulation runs.
# -----------------------------------------------------------------------------

# ----- Data -----

#' Simulate counts from the two-group mixture pi0 f0 + (1 - pi0) f1.
#'
#' A null count is a ZIGP draw. A non-null count is a ZIGP draw plus a shifted
#' draw from f1; with a binomial f1 a non-null count at or below the cutoff is
#' drawn again.
#'
#' @param n Number of counts.
#' @param eta,lambda,theta Parameters of the ZIGP null f0.
#' @param pi0 Proportion of null counts.
#' @param f1 Alternative: "binomial" (Binomial(250, 0.2) + C - 39),
#'   "geometric" (Geometric + C + 1), "binomial mixture" or "poisson mixture".
#' @param C Cutoff: non-null counts lie above it.
#' @param seed Seed, set at the start.
#' @param geometric_prob Success probability of the geometric f1.
#' @return List with `draw` (the counts) and `true_grp` (1 for non-null).
rmixture <- function(n, eta, lambda, theta, pi0, f1 = "binomial", C, seed, geometric_prob = 0.25) {
  set.seed(seed)
  stopifnot(n > 0 || is.na(n) == FALSE)
  stopifnot(eta >= 0 && eta <= 1)
  stopifnot(lambda > 0)
  stopifnot(theta >= 0 && theta <= 1)
  stopifnot(pi0 >= 0 && pi0 <= 1)
  first_component <- 0.6  # weight of the first component of a mixture f1

  draw <- numeric(n)
  true_grp <- numeric(n)
  j <- 0
  while (j < n) {
    if (runif(1) >= pi0) {
      keep <- TRUE
      if (f1 == "binomial") {
        count <- rzigp(eta, lambda, theta) + rbinom(1, 250, 0.20) + C - 39
        keep <- count > C
      } else if (f1 == "geometric") {
        count <- rzigp(eta, lambda, theta) + rgeom(1, geometric_prob) + C + 1
      } else if (f1 == "binomial mixture") {
        if (runif(1) >= first_component) {
          count <- rzigp(eta, lambda, theta) + rbinom(1, 250, 0.2) + C + 1 - 40
          keep <- count > C
        } else {
          count <- rzigp(eta, lambda, theta) + rbinom(1, 50, 0.2) + C + 1 + 0
        }
      } else if (f1 == "poisson mixture") {
        mean <- if (runif(1) >= first_component) 5 else 10
        count <- rzigp(eta, lambda, theta) + rpois(1, mean) + C + 1
      } else {
        stop("unknown f1: ", f1)
      }
      if (keep) {
        draw[j + 1] <- count
        true_grp[j + 1] <- 1
        j <- j + 1
      }
    } else {
      draw[j + 1] <- rzigp(eta, lambda, theta)
      true_grp[j + 1] <- 0
      j <- j + 1
    }
  }
  list(draw = draw, true_grp = true_grp)
}

#' One simulated data set of a setting, split into null and non-null counts.
#'
#' @param setting One row of `simulation_settings()`.
#' @param seed Seed of the data set.
#' @return List: `null` and `alt` (the counts of each group, in the order
#'   drawn), `x` (null counts followed by non-null counts) and `values`,
#'   `null_frequency`, `alt_frequency` (distinct counts and the number of
#'   null and non-null observations at each).
simulated_data <- function(setting, seed) {
  data <- rmixture(setting$N, setting$eta, setting$lambda, setting$theta, setting$pi0,
                   setting$f1, setting$C, seed)
  null <- data$draw[data$true_grp == 0]
  alt <- data$draw[data$true_grp == 1]
  if (length(null) < 2 || length(alt) < 2) stop("a group has fewer than two counts")
  x <- c(null, alt)
  values <- sort(unique(x))
  list(null = null, alt = alt, x = x, values = values,
       null_frequency = as.numeric(table(factor(null, levels = values))),
       alt_frequency = as.numeric(table(factor(alt, levels = values))))
}

# ----- Settings of the simulation runs -----

# The four simulation runs behind the paper. `version` is the version of the
# sampler (SAMPLER_VERSIONS), `bayes_seed` and `eb_seed` the seeds of the data
# sets of the Bayesian and of the empirical Bayes methods (replication d uses
# seed + d), `lambda_init` and `theta_init` starting values of the sampler.
SIMULATION_RUNS <- list(
  run0207 = list(version = "2021-02-07", bayes_seed = 889191, eb_seed = 889191,
                 lambda_init = 2, theta_init = 0.20),
  run0219 = list(version = "2021-02-19", bayes_seed = 889191, eb_seed = 889191,
                 lambda_init = 2, theta_init = 0.25),
  run0224 = list(version = "2021-02-19", bayes_seed = 1440, eb_seed = 889191,
                 lambda_init = 2, theta_init = 0.25),
  run0227 = list(version = "2021-02-19", bayes_seed = 889191, eb_seed = 889191,
                 lambda_init = 2, theta_init = 0.25)
)

# Dirichlet concentrations of the prior sensitivity analysis (Table 2).
PRIOR_GAMMAS <- c(Objective = 0.02, "Jeffreys'" = 0.5, Uniform = 1, "Near Uniform" = 1.5,
                  Informative = 3)

# The six methods, as named in the archive, and the Dirichlet concentration
# (NP, SP) or unused prior argument (GP) of the Bayesian ones in Table 1.
METHODS <- c("NP", "GP", "SP", "TS", "OS", "Storey")
METHOD_LABELS <- c(NP = "Nonparametric", GP = "Parametric", SP = "Semiparametric",
                   TS = "Two-Stage", OS = "One-Stage", Storey = "Storey")
DEFAULT_GAMMA <- c(NP = 1, GP = 0.5, SP = 1)

#' Settings of the simulation runs: one row per run, setting number and N.
#'
#' Every setting has eta = 0.4, pi0 = 0.8 or 0.35, a shifted binomial f1 and
#' 1000 replications. The setting numbers are those of the settings file of
#' each run.
#'
#' @return Data frame with run, setting, eta, lambda, theta, pi0, f1, C, N.
simulation_settings <- function() {
  grid <- function(run, setting, lambda, theta, C, N, pi0 = 0.8) {
    data.frame(run = run, setting = setting, eta = 0.4, lambda = lambda, theta = theta,
               pi0 = pi0, f1 = "binomial", C = C, N = N)
  }
  # run0207: one setting, many N. Settings 16 and 18: Figure S1; 10 and 14: Figure S3.
  n_all <- c(50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
  n_table <- c(50, 100, 200, 500, 1000)
  run0207 <- rbind(grid("run0207", 10, 4.0, 0.30, 5, n_table),
                   grid("run0207", 14, 1.0, 0.20, 5, n_table),
                   grid("run0207", 16, 1.5, 0.00, 10, n_all, pi0 = 0.35),
                   grid("run0207", 18, 1.5, 0.00, 10, n_all))
  # run0219: settings 1 to 30, C = 5; six (lambda, theta) pairs for each N.
  pairs <- data.frame(lambda = c(3, 4, 3, 4, 2, 3), theta = c(0.15, 0.15, 0.20, 0.20, 0.30, 0.30))
  run0219 <- do.call(rbind, lapply(seq_along(c(100, 500, 1000, 50, 200)), function(block) {
    grid("run0219", 6 * (block - 1) + 1:6, pairs$lambda, pairs$theta, 5,
         c(100, 500, 1000, 50, 200)[block])
  }))
  # run0224: settings 1 to 11, C = 10, lambda = 2, theta = 0.3.
  run0224 <- grid("run0224", 1:11, 2, 0.30, 10, n_all)
  # run0227: settings 1 to 20, C = 5; four (lambda, theta) pairs for each N.
  pairs <- data.frame(lambda = c(1, 2, 2, 1), theta = c(0.15, 0.15, 0.20, 0.30))
  run0227 <- do.call(rbind, lapply(seq_along(n_table), function(block) {
    grid("run0227", 4 * (block - 1) + 1:4, pairs$lambda, pairs$theta, 5, n_table[block])
  }))
  rbind(run0207, run0219, run0224, run0227)
}

#' Rows of `simulation_settings()` behind one table or figure.
#'
#' @param item "table1", "table2", "figS1" or "figS3".
#' @return Data frame of settings.
item_settings <- function(item) {
  all <- simulation_settings()
  n_table <- c(50, 100, 200, 500, 1000)
  keep <- switch(item,
    table1 = (all$run == "run0224" & all$N %in% n_table) |
      (all$run == "run0219" & all$lambda == 2 & all$theta == 0.3),
    table2 = all$run == "run0224" & all$N %in% n_table,
    figS1 = all$run == "run0207" & all$setting %in% c(16, 18),
    figS3 = (all$run == "run0207" & all$setting %in% c(10, 14)) | all$run %in% c("run0219", "run0227"),
    stop("unknown item: ", item))
  all[keep, ]
}

# ----- One replication -----

#' Rejections of a Bayesian method in one simulated data set.
#'
#' @param setting One row of `simulation_settings()`.
#' @param model "NP", "SP", "P" or "GP".
#' @param gamma Dirichlet concentration.
#' @param d Replication number; the data set has seed `bayes_seed + d`.
#' @param settings Settings of the paper (`alpha`, `sim_draws`, `sim_burn`).
#' @return One-row data frame: d, R (rejections), V (false rejections), TP
#'   (true rejections), n_alt (non-null counts).
bayes_replication <- function(setting, model, gamma, d, settings) {
  run <- SIMULATION_RUNS[[setting$run]]
  data <- simulated_data(setting, run$bayes_seed + d)
  # No second seed: as in the original runs the sampler continues the stream
  # that simulated_data() has just seeded, so seed + d fixes the data and the chain.
  fit <- bayes_local_fdr(data$x, model, gamma, sampler_settings(run$lambda_init, run$theta_init),
                         draws = settings$sim_draws, burn = settings$sim_burn,
                         version = run$version)
  rejected <- which(fit$fdr <= settings$alpha)
  V <- sum(data$null_frequency[rejected])
  TP <- sum(data$alt_frequency[rejected])
  data.frame(d = d, R = TP + V, V = V, TP = TP, n_alt = length(data$alt))
}

#' Rejections of the three empirical Bayes methods in one simulated data set.
#'
#' @param setting One row of `simulation_settings()`.
#' @param d Replication number; the data set has seed `eb_seed + d`.
#' @param settings Settings of the paper (`alpha`).
#' @return One-row data frame: d, n_alt, and for each of TS (two-stage), OS
#'   (one-stage) and Storey the rejections `_R`, false rejections `_V` and
#'   true rejections `_TP`.
eb_replication <- function(setting, d, settings) {
  data <- simulated_data(setting, SIMULATION_RUNS[[setting$run]]$eb_seed + d)
  data.frame(d = d, n_alt = length(data$alt),
             eb_simulation_counts(data, setting$C, settings$alpha))
}

# ----- Driver -----

#' Folder of the stored replications of a setting and method.
#'
#' Settings of different runs that share data sets and sampler version share
#' their replications.
#'
#' @param setting One row of `simulation_settings()`.
#' @param method "EB", or a Bayesian model.
#' @param gamma Dirichlet concentration (Bayesian methods).
#' @param settings Settings of the paper.
#' @param results Results folder.
#' @return Path of the folder.
replication_dir <- function(setting, method, gamma, settings, results) {
  run <- SIMULATION_RUNS[[setting$run]]
  data <- sprintf("eta%g_lambda%g_theta%g_pi0%g_%s_C%g_N%g", setting$eta, setting$lambda,
                  setting$theta, setting$pi0, setting$f1, setting$C, setting$N)
  how <- if (method == "EB") {
    sprintf("EB_seed%d", run$eb_seed)
  } else {
    sprintf("%s_gamma%g_seed%d_sampler%s_init%g_draws%d_burn%d", method, gamma, run$bayes_seed,
            run$version, run$theta_init, settings$sim_draws, settings$sim_burn)
  }
  file.path(results, "replications", data, how)
}

#' Run the replications of a set of settings that are not yet stored.
#'
#' One csv file is written per setting, method and replication, so a study
#' can be run in pieces and continued.
#'
#' @param cells Data frame of settings with the columns `method` ("EB" or a
#'   Bayesian model) and `gamma`.
#' @param runs Replication numbers to run.
#' @param settings Settings of the paper.
#' @param results Results folder.
#' @param cores Number of parallel processes.
#' @return NULL, invisibly.
run_replications <- function(cells, runs, settings, results, cores = 1) {
  jobs <- expand.grid(cell = seq_len(nrow(cells)), d = runs)
  jobs$file <- file.path(vapply(jobs$cell, function(i) {
    replication_dir(cells[i, ], cells$method[i], cells$gamma[i], settings, results)
  }, character(1)), sprintf("rep%04d.csv", jobs$d))
  jobs <- jobs[!file.exists(jobs$file), ]
  log_message("%d replications to run, %d cores", nrow(jobs), cores)
  status <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
    cell <- cells[jobs$cell[i], ]
    result <- if (cell$method == "EB") {
      eb_replication(cell, jobs$d[i], settings)
    } else {
      bayes_replication(cell, cell$method, cell$gamma, jobs$d[i], settings)
    }
    dir.create(dirname(jobs$file[i]), showWarnings = FALSE, recursive = TRUE)
    write.csv(result, jobs$file[i], row.names = FALSE)
    TRUE
  }, mc.cores = cores, mc.preschedule = FALSE)
  failed <- !vapply(status, isTRUE, logical(1))
  if (any(failed)) stop(sum(failed), " replications failed; first: ", status[[which(failed)[1]]])
  invisible(NULL)
}

#' Mean and standard deviation of R, FDR and TPR over stored replications.
#'
#' Rounded as in the paper: R to two decimals, FDR and TPR to four. The FDR of
#' a replication without rejections is 0 (Storey: left out of the mean).
#'
#' @param cells,settings,results As in `run_replications()`.
#' @param replications Number of replications D; replications 1 to D are read.
#' @return Data frame with run, setting, N, method, gamma and R, R_sd, FDR,
#'   FDR_sd, TPR, TPR_sd; one row per Bayesian cell and three per EB cell.
summarize_replications <- function(cells, replications, settings, results) {
  summary_row <- function(cell, method, gamma, R, fdr, tpr) {
    data.frame(run = cell$run, setting = cell$setting, N = cell$N, method = method, gamma = gamma,
               R = round(mean(R, na.rm = TRUE), 2), R_sd = round(sd(R, na.rm = TRUE), 2),
               FDR = round(mean(fdr, na.rm = TRUE), 4), FDR_sd = round(sd(fdr, na.rm = TRUE), 4),
               TPR = round(mean(tpr, na.rm = TRUE), 4), TPR_sd = round(sd(tpr, na.rm = TRUE), 4))
  }
  do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    cell <- cells[i, ]
    files <- file.path(replication_dir(cell, cell$method, cell$gamma, settings, results),
                       sprintf("rep%04d.csv", seq_len(replications)))
    if (!all(file.exists(files))) stop("replications missing for ", dirname(files[1]))
    reps <- do.call(rbind, lapply(files, read.csv))
    zero_if_none <- function(V, R) ifelse(R == 0, 0, V / R)
    if (cell$method != "EB") {
      return(summary_row(cell, cell$method, cell$gamma, reps$R, zero_if_none(reps$V, reps$R),
                         reps$TP / reps$n_alt))
    }
    rbind(summary_row(cell, "TS", NA, reps$TS_R, zero_if_none(reps$TS_V, reps$TS_R),
                      reps$TS_TP / reps$n_alt),
          summary_row(cell, "OS", NA, reps$OS_R, zero_if_none(reps$OS_V, reps$OS_R),
                      reps$OS_TP / reps$n_alt),
          summary_row(cell, "Storey", NA, reps$Storey_R, reps$Storey_V / reps$Storey_R,
                      reps$Storey_TP / reps$n_alt))
  }))
}

#' Summaries of the settings of a table or figure, by simulation.
#'
#' Runs the replications asked for that are not stored yet; when replications
#' 1 to `opts$replications` are all stored, returns their summaries.
#'
#' @param item "table1", "table2", "figS1" or "figS3".
#' @param opts Options: `runs` (replications to run now), `replications` (D),
#'   `cores`.
#' @param settings Settings of the paper.
#' @param paths Project paths.
#' @return Data frame as `summarize_replications()`, or NULL if replications
#'   are still missing.
simulation_summaries <- function(item, opts, settings, paths) {
  setup <- item_settings(item)
  with_method <- function(method, gamma) cbind(setup, method = method, gamma = gamma)
  cells <- if (item == "table2") {
    do.call(rbind, lapply(PRIOR_GAMMAS, function(g) rbind(with_method("NP", g), with_method("SP", g))))
  } else {
    rbind(with_method("NP", DEFAULT_GAMMA[["NP"]]), with_method("GP", DEFAULT_GAMMA[["GP"]]),
          with_method("SP", DEFAULT_GAMMA[["SP"]]), with_method("EB", NA))
  }
  runs <- parse_index_set(opts$runs, count = TRUE)
  run_replications(cells, runs[runs <= opts$replications], settings, paths$results, opts$cores)
  done <- vapply(seq_len(nrow(cells)), function(i) {
    all(file.exists(file.path(replication_dir(cells[i, ], cells$method[i], cells$gamma[i],
                                              settings, paths$results),
                              sprintf("rep%04d.csv", seq_len(opts$replications)))))
  }, logical(1))
  if (!all(done)) {
    log_message("replications 1 to %d are not all stored yet; no summary written", opts$replications)
    return(NULL)
  }
  summarize_replications(cells, opts$replications, settings, paths$results)
}

#' Archived summaries of the original runs for a table or figure.
#'
#' @param item "table1", "table2", "figS1" or "figS3".
#' @param paths Project paths.
#' @return Data frame with run, setting, N, method, gamma and the archived
#'   R, R_sd, FDR, FDR_sd, TPR, TPR_sd (NA where a run did not store them).
archived_summaries <- function(item, paths) {
  archive <- read.csv(file.path(paths$archive, "simulation_summaries.csv"))
  setup <- item_settings(item)
  archive <- archive[paste(archive$run, archive$setting, archive$N) %in%
                       paste(setup$run, setup$setting, setup$N), ]
  if (item == "table2") archive[archive$study == "prior", ] else archive[archive$study == "methods", ]
}

#' Add the parameters of the settings to a table of summaries.
#'
#' @param summaries Data frame with run, setting and N.
#' @return `summaries` with the columns lambda, theta, pi0 and C.
with_parameters <- function(summaries) {
  all <- simulation_settings()
  row <- match(paste(summaries$run, summaries$setting, summaries$N), paste(all$run, all$setting, all$N))
  cbind(summaries, all[row, c("lambda", "theta", "pi0", "C")])
}

# ----- Layout of the tables -----

#' The six summary columns of a block of a table, with a prefix.
#' @param summaries Rows of summaries, in the order of the table rows.
#' @param prefix Prefix of the column names.
#' @return Data frame with <prefix>_R, _R_sd, _FDR, _FDR_sd, _TPR, _TPR_sd.
summary_block <- function(summaries, prefix) {
  block <- summaries[c("R", "R_sd", "FDR", "FDR_sd", "TPR", "TPR_sd")]
  setNames(block, paste(prefix, names(block), sep = "_"))
}

#' Table 1: R, FDR and TPR of the six methods by N, for C = 10 and C = 5.
#' @param summaries Summaries of the settings of `item_settings("table1")`.
#' @return Data frame with one row per N and method.
table1_layout <- function(summaries) {
  summaries <- with_parameters(summaries)
  rows <- expand.grid(method = METHODS, N = sort(unique(summaries$N)), stringsAsFactors = FALSE)
  block <- function(C) {
    summaries[match(paste(rows$N, rows$method, C), paste(summaries$N, summaries$method, summaries$C)), ]
  }
  data.frame(N = rows$N, method = METHOD_LABELS[rows$method], summary_block(block(10), "C10"),
             summary_block(block(5), "C5"), row.names = NULL)
}

#' Table 2: R, FDR and TPR of NP and SP by N and prior.
#' @param summaries Summaries of the settings of `item_settings("table2")`.
#' @return Data frame with one row per N and prior.
table2_layout <- function(summaries) {
  rows <- expand.grid(gamma = PRIOR_GAMMAS, N = sort(unique(summaries$N)))
  block <- function(method) {
    summaries[match(paste(rows$N, rows$gamma, method),
                    paste(summaries$N, summaries$gamma, summaries$method)), ]
  }
  data.frame(N = rows$N, prior = names(PRIOR_GAMMAS)[match(rows$gamma, PRIOR_GAMMAS)],
             gamma = rows$gamma, summary_block(block("NP"), "NP"), summary_block(block("SP"), "SP"),
             row.names = NULL)
}
