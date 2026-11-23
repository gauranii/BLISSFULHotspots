# -----------------------------------------------------------------------------
# File:    checks/benchmark_chain.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Time one Gibbs chain at the settings of Table 3 (10000 stored
#          iterations after a burn-in of 1000, sampler version 2018-03-25,
#          seed 1234), in R or with the C++ functions of rcpp/src/. Only the
#          chain is timed: not the start of R, not the loading of the code.
# Inputs:  data/domains/<domain>.txt; rcpp/src/ with rcpp.
# Outputs: One csv line on standard output: impl,domain,model,seconds.
#          The local fdr of every count of the chain to <out file>, so that
#          the R and the Rcpp chains can be compared.
# Usage:   Rscript checks/benchmark_chain.R r|rcpp <domain> <model> <out file>
#          Run from the repository root, inside `nix develop`.
# Source:  domain_bayes_rejections() in R/real_data.R.
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
impl <- args[1]
domain <- args[2]
model <- args[3]
out <- args[4]

source("R/setup.R")
if (impl == "rcpp") source("rcpp/R/load.R")

x <- read_domain(domain, PATHS$domains)$counts
seconds <- system.time({
  set.seed(PAPER$real_seed)
  fit <- bayes_local_fdr(x, model, PAPER$real_gamma[[model]],
                         sampler_settings(lambda_init = 1, theta_init = 0.5),
                         draws = PAPER$real_draws, burn = PAPER$real_burn,
                         version = "2018-03-25")
})[["elapsed"]]
write.csv(data.frame(count = fit$values, fdr = sprintf("%.17g", fit$fdr)), out, row.names = FALSE)
cat(sprintf("%s,%s,%s,%.3f\n", impl, domain, model, seconds))
