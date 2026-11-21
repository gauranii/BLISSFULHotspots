// -----------------------------------------------------------------------------
// File:    rcpp/src/distributions.cpp
// Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
//          false discovery rate for sparse count data with application to the
//          discovery of hotspots in protein domains. Ann. Appl. Stat.
//          16(3):1459-1475. Tables 1 to 3, Figures S1 to S3.
// Purpose: The functions of R/distributions.R as R functions written in C++:
//          same names, same arguments, same values and the same random
//          numbers drawn. The code itself is in distributions_core.h.
// Inputs:  None.
// Outputs: None.
// Usage:   source("R/setup.R"); source("rcpp/R/load.R")
//          dzigp(0:5, eta = 0.4, lambda = 2, theta = 0.3)
// Source:  R/distributions.R.
// -----------------------------------------------------------------------------
// [[Rcpp::plugins(localfdr)]]
#include "distributions_core.h"

using localfdr::dvec;

// [[Rcpp::export]]
Rcpp::NumericVector dzigp(Rcpp::NumericVector x, double eta, double lambda, double theta,
                          bool check = true) {
  return Rcpp::wrap(localfdr::dzigp(Rcpp::as<dvec>(x), eta, lambda, theta, check));
}

// [[Rcpp::export]]
double rzigp(double eta, double lambda, double theta, bool check = true) {
  return localfdr::rzigp(eta, lambda, theta, check);
}

// [[Rcpp::export]]
Rcpp::NumericVector rdirichlet(Rcpp::NumericVector beta, bool replace_missing = false) {
  return Rcpp::wrap(localfdr::rdirichlet(Rcpp::as<dvec>(beta), replace_missing));
}

// [[Rcpp::export]]
Rcpp::NumericVector rmultinomial(double size, Rcpp::NumericVector prob) {
  return Rcpp::wrap(localfdr::rmultinomial(size, Rcpp::as<dvec>(prob)));
}

// [[Rcpp::export]]
Rcpp::NumericVector shifted_density(std::string model, Rcpp::NumericVector x, double xi,
                                    double delta, double nu, double C) {
  return Rcpp::wrap(localfdr::shifted_density(localfdr::model_code(model), Rcpp::as<dvec>(x),
                                              xi, delta, nu, C));
}
