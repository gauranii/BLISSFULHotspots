# -----------------------------------------------------------------------------
# File:    R/real_data.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Table 3, Figure 1.
# Purpose: Reader of the protein domain files and the analysis of a domain:
#          the number of positions declared mutational hotspots by the
#          empirical Bayes method and by the Bayesian methods.
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R")
# Source:  real_scenarios_wz.R, real_scenarios.R (gibbs_out, RealDataFDR),
#          fcns6_initial.R; protein.R (reading of the files for Figure 1).
# -----------------------------------------------------------------------------

# The six protein domains, in the order of Figure 1.
DOMAINS <- c("cd00054", "cd00083", "cd00590", "smart00179", "smart00209", "cd00063")

#' Read the mutation counts of a protein domain.
#'
#' A domain file has the number of positions N on its first line, then one
#' line `position:count:mutations` per position, starting with a position 0
#' that holds no mutation.
#'
#' @param domain Name of the domain, for example "cd00083".
#' @param folder Folder of the domain files.
#' @return List with `N`, `counts` (mutation counts of positions 1 to N) and
#'   `with_position_0` (the same with the count of position 0 in front).
read_domain <- function(domain, folder) {
  table <- read.delim(file.path(folder, paste0(domain, ".txt")), sep = ":", header = FALSE)
  list(N = table[1, 1], counts = table[3:nrow(table), 2], with_position_0 = table[2:nrow(table), 2])
}

#' Number of positions at each mutation count 0, 1, ..., max.
#'
#' @param x Mutation counts.
#' @return Named vector of frequencies.
count_frequencies <- function(x) {
  frequency <- c(sum(x == 0), tabulate(x, max(x)))
  names(frequency) <- 0:max(x)
  frequency
}

#' Hotspot positions of a domain by one Bayesian method.
#'
#' Runs the sampler as in the analysis of the paper (seed, iterations and
#' Dirichlet concentration from `settings`) and counts the positions whose
#' mutation count has posterior mean local fdr at most alpha.
#'
#' @param x Mutation counts of the positions.
#' @param model "NP", "SP", "P" or "GP".
#' @param settings Settings of the paper.
#' @return Number of positions declared hotspots.
domain_bayes_rejections <- function(x, model, settings) {
  set.seed(settings$real_seed)
  fit <- bayes_local_fdr(x, model, settings$real_gamma[[model]],
                         sampler_settings(lambda_init = 1, theta_init = 0.5),
                         draws = settings$real_draws, burn = settings$real_burn,
                         version = "2018-03-25")
  sum(fit$frequency[which(fit$fdr <= settings$alpha)])
}

#' Hotspot positions of the domains by every method (Table 3).
#'
#' @param domains Names of the domains.
#' @param models Bayesian models to run.
#' @param settings Settings of the paper.
#' @param paths Project paths.
#' @param cores Number of parallel processes.
#' @return Data frame with domain, N, EB and one column per model.
domain_rejections <- function(domains, models, settings, paths, cores = 1) {
  jobs <- expand.grid(model = models, domain = domains, stringsAsFactors = FALSE)
  counts <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
    x <- read_domain(jobs$domain[i], paths$domains)$counts
    domain_bayes_rejections(x, jobs$model[i], settings)
  }, mc.cores = cores, mc.preschedule = FALSE)
  failed <- !vapply(counts, is.numeric, logical(1))
  if (any(failed)) stop("a chain failed: ", counts[[which(failed)[1]]])
  bayes <- matrix(unlist(counts), ncol = length(models), byrow = TRUE, dimnames = list(NULL, models))
  eb <- vapply(domains, function(d) eb_rejections(read_domain(d, paths$domains)$counts, settings$alpha),
               numeric(1))
  N <- vapply(domains, function(d) length(read_domain(d, paths$domains)$counts), numeric(1))
  data.frame(domain = domains, N = N, EB = eb, bayes, row.names = NULL)
}
