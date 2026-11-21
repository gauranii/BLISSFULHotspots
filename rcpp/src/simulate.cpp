// -----------------------------------------------------------------------------
// File:    rcpp/src/simulate.cpp
// Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
//          false discovery rate for sparse count data with application to the
//          discovery of hotspots in protein domains. Ann. Appl. Stat.
//          16(3):1459-1475. Tables 1 and 2, Figures S1 to S3.
// Purpose: rmixture() of R/simulate.R in C++: simulated counts from the
//          two-group mixture pi0 f0 + (1 - pi0) f1, with the same name and
//          arguments, the same seed set at the start and the same random
//          numbers drawn in the same order. The rest of R/simulate.R (the
//          settings of the runs and the driver of the replications) stays
//          in R and calls this function.
// Inputs:  None.
// Outputs: None.
// Usage:   source("R/setup.R"); source("rcpp/R/load.R")
//          rmixture(100, eta = 0.4, lambda = 2, theta = 0.3, pi0 = 0.8,
//                   f1 = "binomial", C = 5, seed = 1)
// Source:  R/simulate.R (rmixture).
// -----------------------------------------------------------------------------
// [[Rcpp::plugins(localfdr)]]
#include "distributions_core.h"

// [[Rcpp::export]]
Rcpp::List rmixture(double n, double eta, double lambda, double theta, double pi0,
                    std::string f1, double C, SEXP seed, double geometric_prob = 0.25) {
  Rcpp::Function set_seed("set.seed", R_BaseEnv);
  set_seed(seed);
  if (ISNAN(n)) Rcpp::stop("n > 0 || is.na(n) == FALSE is not TRUE");
  if (n < 0) Rcpp::stop("invalid 'length' argument");
  localfdr::check_zigp(eta, lambda, theta);
  if (!(pi0 >= 0 && pi0 <= 1)) Rcpp::stop("pi0 >= 0 && pi0 <= 1 is not TRUE");
  const double first_component = 0.6;  // weight of the first component of a mixture f1

  std::size_t size = (std::size_t) n;
  Rcpp::NumericVector draw(size), true_grp(size);
  std::size_t j = 0;
  while ((double) j < n) {
    if (R::runif(0.0, 1.0) >= pi0) {
      bool keep = true;
      double count;
      if (f1 == "binomial") {
        count = localfdr::rzigp(eta, lambda, theta, true);
        count = count + localfdr::rbinom_one(250, 0.20) + C - 39;
        keep = count > C;
      } else if (f1 == "geometric") {
        count = localfdr::rzigp(eta, lambda, theta, true);
        count = count + R::rgeom(geometric_prob) + C + 1;
      } else if (f1 == "binomial mixture") {
        if (R::runif(0.0, 1.0) >= first_component) {
          count = localfdr::rzigp(eta, lambda, theta, true);
          count = count + localfdr::rbinom_one(250, 0.2) + C + 1 - 40;
          keep = count > C;
        } else {
          count = localfdr::rzigp(eta, lambda, theta, true);
          count = count + localfdr::rbinom_one(50, 0.2) + C + 1 + 0;
        }
      } else if (f1 == "poisson mixture") {
        double mean = R::runif(0.0, 1.0) >= first_component ? 5 : 10;
        count = localfdr::rzigp(eta, lambda, theta, true);
        count = count + R::rpois(mean) + C + 1;
      } else {
        Rcpp::stop("unknown f1: %s", f1);
      }
      if (ISNAN(count)) localfdr::stop_missing_condition();
      if (keep) {
        if (j >= size) Rcpp::stop("subscript out of bounds");
        draw[j] = count;
        true_grp[j] = 1;
        j = j + 1;
      }
    } else {
      if (j >= size) Rcpp::stop("subscript out of bounds");
      draw[j] = localfdr::rzigp(eta, lambda, theta, true);
      true_grp[j] = 0;
      j = j + 1;
    }
  }
  return Rcpp::List::create(Rcpp::Named("draw") = draw, Rcpp::Named("true_grp") = true_grp);
}
