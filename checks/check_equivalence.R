# -----------------------------------------------------------------------------
# File:    checks/check_equivalence.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Run the original functions (sourced unchanged from the folders of
#          the original runs) and the functions of R/ on the same small cases
#          and report whether their results are identical: simulated data,
#          two replications of a Table 1 setting for every method, one
#          replication of the other simulation runs and of two priors of
#          Table 2, shortened chains of one protein domain for the four
#          Bayesian models with its empirical Bayes count, and one panel of
#          Figure S2.
# Inputs:  The run folders in the folder named by R012_RUNS (see
#          archive/reduce_archive.R); data/domains/cd00083.txt.
# Outputs: prints one line per comparison and a summary; exit status 1 unless
#          every comparison is identical.
# Usage:   R012_RUNS=/path/to/run/folders Rscript checks/check_equivalence.R
#          [--cores=2] [--large]
#          --large: also one replication of a setting with N = 1000 for every
#          method. About 1 minute on 2 cores, 2 minutes with --large. The
#          original functions call library() for packages they do not use in
#          these cases; that call is skipped.
# -----------------------------------------------------------------------------

source("R/setup.R")

opts <- parse_args(list(cores = 2, large = FALSE))
runs_dir <- Sys.getenv("R012_RUNS")
if (!nzchar(runs_dir))
  stop("Set R012_RUNS to the folder of the original runs (not part of this repository).")
stopifnot(dir.exists(runs_dir))

#' Source a file of an original run, and the files it sources, into a new
#' environment.
#' @param ... Path of the file under the run folders.
#' @return The environment.
original <- function(...) {
  env <- new.env()
  env$source <- function(file, ...) base::source(file, local = env)
  env$library <- function(...) invisible(NULL)
  file <- file.path(runs_dir, ...)
  old <- setwd(dirname(file))
  on.exit(setwd(old))
  base::source(basename(file), local = env)
  env
}
ORIGINAL_SETTINGS <- c(run0207 = "2021-02-07 Figures S1 and S3 simulations (RERUN 0216)/CODES/initial_vals0207.R",
                       run0219 = "2021-02-19 Table 1 (C = 5) and Figure S3 simulations (EXCEL 0220)/sim_initial0219.R",
                       run0224 = "2021-02-24 Table 1 (C = 10) simulations (RERUN 0225)/sim_initial0225.R",
                       run0227 = "2021-02-25 Table 2 and Figure S3 simulations (RERUN 0309)/sim_initial0227.R")
ALGORITHM <- c(P = 1, GP = 2, SP = 3, NP = 4)
quietly <- function(expr) {
  sink(nullfile())
  on.exit(sink())
  suppressWarnings(expr)
}
setting_of <- function(run, setting, N) {
  all <- simulation_settings()
  all[all$run == run & all$setting == setting & all$N == N, ]
}

#' Compare one simulation setting: replications 1 to D of a method.
#' @return Data frame with one row per comparison: `case` and `identical`.
simulation_case <- function(run, setting, N, method, gamma, D) {
  cell <- setting_of(run, setting, N)
  env <- original(ORIGINAL_SETTINGS[[run]])
  truth <- list(n = N, D = D, eta_true = cell$eta, lambda_true = cell$lambda,
                theta_true = cell$theta, pi0_true = cell$pi0, B = 1, C_true = cell$C)
  if (method == "EB") {
    old <- quietly(do.call(env$EM_fdr, c(truth, list(FDRLevel = 0.05))))
    new <- do.call(rbind, lapply(seq_len(D), function(d) eb_replication(cell, d, PAPER)))
    zero_if_none <- function(V, R) ifelse(R == 0, 0, V / R)
    same <- c(
      "two-stage R, FDR, TPR" = identical(old$R[1, 1, ], as.numeric(new$TS_R)) &&
        identical(old$FDR[1, 1, ], zero_if_none(new$TS_V, new$TS_R)) &&
        identical(old$TPR[1, 1, ], new$TS_TP / new$n_alt),
      "one-stage R, FDR, TPR" = identical(old$R[2, 1, ], as.numeric(new$OS_R)) &&
        identical(old$FDR[2, 1, ], zero_if_none(new$OS_V, new$OS_R)) &&
        identical(old$TPR[2, 1, ], new$OS_TP / new$n_alt),
      # The original returns Storey's procedure only as rounded means.
      "Storey mean R, FDR, TPR" = identical(
        as.numeric(old$summary[1, 7:9]),
        c(round(mean(new$Storey_R), 2), round(mean(new$Storey_V / new$Storey_R, na.rm = TRUE), 4),
          round(mean(new$Storey_TP / new$n_alt), 4))))
  } else {
    old <- quietly(do.call(env$bayes_fdr, c(truth, list(
      alpha = env$alpha, zigp = TRUE, drawC = TRUE, algorithm = ALGORITHM[[method]], gamma = gamma,
      phi_init = env$phi_init, hyper_tau = env$hyper_tau, hyper_xi = env$hyper_xi,
      proposal = env$proposal, draw = 5000, thin = 1, burn = 1000, report_period = 6001))))
    new <- do.call(rbind, lapply(seq_len(D), function(d) bayes_replication(cell, method, gamma, d, PAPER)))
    same <- c("R, FDR, TPR" = identical(old$R, as.numeric(new$R)) &&
                identical(old$FDR, ifelse(new$R == 0, 0, new$V / new$R)) &&
                identical(old$TPR, new$TP / new$n_alt))
  }
  label <- sprintf("%s setting %d (N = %d), %s%s, replications 1 to %d: %s", run, setting, N, method,
                   if (method == "EB") "" else sprintf(" gamma %g", gamma), D, names(same))
  data.frame(case = label, identical = unname(same))
}

#' Compare shortened chains of one protein domain, and its EB count.
domain_case <- function(model) {
  folder <- "2018-03-25 Table 3 protein domain analysis"
  x <- read_domain("cd00083", PATHS$domains)$counts
  if (model == "EB") {
    env <- original(folder, "real_scenarios.R")
    old <- quietly(env$RealDataFDR(x, 0.05)$R[1, 1])
    return(data.frame(case = "cd00083, EB: positions declared",
                      identical = identical(as.numeric(old), as.numeric(eb_rejections(x, 0.05)))))
  }
  env <- original(folder, "real_scenarios.R")
  gamma <- PAPER$real_gamma[[model]]
  set.seed(PAPER$real_seed)
  old <- quietly(env$gibbs(model = ALGORITHM[[model]], TRUE, TRUE, x, gamma = gamma, env$phi_init,
                           env$pi0_prior, env$hyper_tau, env$hyper_xi, env$proposal, env$adapt,
                           draw = 2000, thin = 1, burn = 300, report_period = 10000))
  set.seed(PAPER$real_seed)
  new <- bayes_local_fdr(x, model, gamma, sampler_settings(lambda_init = 1, theta_init = 0.5),
                         draws = 2000, burn = 300, version = "2018-03-25")
  data.frame(case = sprintf("cd00083, %s gamma %g, 2300 iterations: fdr and draws of C, pi0, eta, lambda, theta",
                            model, gamma),
             identical = identical(old$fdr, new$fdr) && identical(old$C, new$C) &&
               identical(old$pi0, new$pi0) && identical(old$eta, new$eta) &&
               identical(old$lambda, new$lambda) && identical(old$theta, new$theta))
}

#' Compare simulated data sets and the frequencies of one panel of Figure S2.
data_case <- function() {
  env <- original(ORIGINAL_SETTINGS[["run0224"]])
  f1 <- c("binomial", "geometric", "binomial mixture", "poisson mixture")
  same_data <- vapply(1:4, function(B) {
    identical(env$rmixture(200, 0.4, 2, 0.3, 0.8, B, 10, 1440 + B),
              rmixture(200, 0.4, 2, 0.3, 0.8, f1[B], 10, 1440 + B))
  }, logical(1))
  # Figure S2, N = 200, lambda = 2: the original script as it is.
  folder <- file.path(runs_dir, "2021-05-04 Figure S2 histograms (RERUN 0204)")
  env <- original("2021-05-04 Figure S2 histograms (RERUN 0204)", "simFDR_fcns.R")
  work <- tempfile("figS2")
  dir.create(work)
  old <- setwd(work)
  on.exit(setwd(old))
  quietly(base::source(file.path(folder, "barplots0504.R"), local = env))
  setting <- data.frame(N = 200, eta = 0.4, lambda = 2, theta = 0.3, pi0 = 0.8, f1 = "binomial", C = 5)
  data <- simulated_data(setting, PAPER$figS2_seed)
  new <- mixture_barplot(data$null, data$alt, file.path(work, "new.pdf"), PAPER$figure_size)
  data.frame(case = c(sprintf("simulated data, f1 %s, N = 200", f1),
                      "Figure S2, N = 200, lambda = 2: plotted frequencies"),
             identical = c(same_data, identical(env$counts1, new$null) && identical(env$counts2, new$alt)))
}

cases <- c(
  list(function() data_case()),
  lapply(c("NP", "GP", "SP", "EB"), function(m) {
    function() simulation_case("run0224", 1, 50, m, DEFAULT_GAMMA[m], D = 2)
  }),
  lapply(c("NP", "GP", "SP", "EB"), function(m) {
    function() simulation_case("run0219", 23, 50, m, DEFAULT_GAMMA[m], D = 1)
  }),
  lapply(c("NP", "GP", "SP", "EB"), function(m) {
    function() simulation_case("run0207", 16, 50, m, DEFAULT_GAMMA[m], D = 1)
  }),
  list(function() simulation_case("run0227", 1, 50, "NP", 1, D = 1),
       function() simulation_case("run0224", 1, 50, "NP", 0.02, D = 1),
       function() simulation_case("run0224", 1, 50, "SP", 3, D = 1)),
  lapply(c("EB", "NP", "SP", "P", "GP"), function(m) function() domain_case(m)),
  if (opts$large) lapply(c("NP", "GP", "SP", "EB"), function(m) {
    function() simulation_case("run0224", 11, 1000, m, DEFAULT_GAMMA[m], D = 1)
  })
)
results <- parallel::mclapply(cases, function(case) case(), mc.cores = opts$cores, mc.preschedule = FALSE)
failed <- !vapply(results, is.data.frame, logical(1))
if (any(failed)) stop("a case failed: ", results[[which(failed)[1]]])
results <- do.call(rbind, results)
cat(sprintf("%-9s  %s\n", ifelse(results$identical, "identical", "DIFFERENT"), results$case), sep = "")
cat(sprintf("%d of %d comparisons identical\n", sum(results$identical), nrow(results)))
quit(status = if (all(results$identical)) 0 else 1)
