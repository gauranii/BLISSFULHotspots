// -----------------------------------------------------------------------------
// File:    rcpp/src/distributions_core.h
// Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
//          false discovery rate for sparse count data with application to the
//          discovery of hotspots in protein domains. Ann. Appl. Stat.
//          16(3):1459-1475. Tables 1 to 3, Figures S1 to S3.
// Purpose: The functions of R/distributions.R in C++, shared by the files of
//          rcpp/src/: the zero-inflated generalized Poisson (ZIGP) null
//          distribution, the shifted Poisson and generalized Poisson
//          alternatives, and the Dirichlet and multinomial samplers. Each
//          follows its R function operation by operation and draws from R's
//          random number generator in the same order.
// Inputs:  None.
// Outputs: None.
// Usage:   #include "distributions_core.h"
// Source:  R/distributions.R (dzigp, rzigp, rdirichlet, rmultinomial,
//          shifted_density); R 4.6.1 src/library/stats/src/random.c (how
//          runif(), rgamma() and rbinom() call the generators).
// -----------------------------------------------------------------------------
#ifndef LOCALFDR_DISTRIBUTIONS_CORE_H
#define LOCALFDR_DISTRIBUTIONS_CORE_H

#include "rbase.h"

#include <climits>

namespace localfdr {

// Alternatives f1; OTHER is any other text, which R treats as "no f1".
enum Model { MODEL_P, MODEL_GP, MODEL_SP, MODEL_NP, MODEL_OTHER };

inline Model model_code(const std::string& model) {
  if (model == "P") return MODEL_P;
  if (model == "GP") return MODEL_GP;
  if (model == "SP") return MODEL_SP;
  if (model == "NP") return MODEL_NP;
  return MODEL_OTHER;
}

// The range checks of dzigp() and rzigp().
inline void check_zigp(double eta, double lambda, double theta) {
  if (!(eta >= 0 && eta <= 1)) Rcpp::stop("eta >= 0 && eta <= 1 is not TRUE");
  if (!(lambda > 0)) Rcpp::stop("lambda > 0 is not TRUE");
  if (!(theta >= 0 && theta <= 1)) Rcpp::stop("theta >= 0 && theta <= 1 is not TRUE");
}

// Probability mass function of the ZIGP at the counts x.
inline dvec dzigp(const dvec& x, double eta, double lambda, double theta, bool check) {
  if (check) check_zigp(eta, lambda, theta);
  dvec density(x.size(), 0.0);
  for (std::size_t i = 0; i < x.size(); i++) {
    double gp = r_exp(r_log(lambda) + (x[i] - 1) * r_log(lambda + theta * x[i]) -
                      r_lgamma(x[i] + 1) - (lambda + theta * x[i]));
    if (ISNAN(x[i])) continue;  // which() leaves a missing count out: its density stays 0
    density[i] = x[i] == 0 ? eta + (1 - eta) * gp : (1 - eta) * gp;
  }
  return density;
}

// One draw from the ZIGP by inversion; one uniform number.
inline double rzigp(double eta, double lambda, double theta, bool check) {
  if (check) check_zigp(eta, lambda, theta);
  dvec cdf;
  double unif = R::runif(0.0, 1.0);
  cdf.push_back(eta + (1 - eta) * r_exp(-lambda));

  if (unif > 0) {
    bool recursive = false;
    double j = 1, p = 0, log_p_previous = 0;
    for (;;) {
      double current = cdf[(std::size_t) j - 1];
      if (ISNAN(current)) stop_missing_condition();
      if (!(current < unif)) break;
      if (!recursive) {
        p = r_exp(r_log(lambda) + (j - 1) * r_log(lambda + theta * j) - r_lgamma(j + 1) -
                  (lambda + theta * j) + r_log(1 - eta));
      }
      if (ISNAN(p)) stop_missing_condition();
      if (p == R_PosInf) {
        recursive = true;
        log_p_previous = r_log(lambda) + (j - 2) * r_log(lambda + theta * (j - 1)) -
          r_lgamma(j) - (lambda + theta * (j - 1)) + r_log(1 - eta);
      }
      if (recursive) {
        double log_p = log_p_previous + r_log(lambda + theta * (j - 1)) - r_log(j) - theta +
          (j - 1) * r_log(lambda + theta * j) - (j - 1) * r_log(lambda + theta * (j - 1));
        log_p_previous = log_p;
        p = r_exp(log_p);
      }
      cdf.push_back(current + p);
      j = j + 1;
    }
  }

  std::size_t k = 0;
  for (;;) {
    if (k >= cdf.size() || ISNAN(cdf[k])) stop_missing_condition();
    if (!(unif > cdf[k])) break;
    k++;
  }
  return (double) k;
}

// One draw from a Dirichlet distribution, from independent gamma draws.
inline dvec rdirichlet(dvec beta, bool replace_missing) {
  if (replace_missing) {
    for (std::size_t i = 0; i < beta.size(); i++) {
      if (ISNAN(beta[i])) beta[i] = 0.00001;
    }
  }
  double scale = 1 / 1.0;  // rgamma(n, shape, rate = 1) draws with scale 1 / rate
  dvec draw(beta.size());
  for (std::size_t i = 0; i < beta.size(); i++) draw[i] = R::rgamma(beta[i], scale);
  double total = 0;
  for (std::size_t i = 0; i < draw.size(); i++) total = total + draw[i];
  for (std::size_t i = 0; i < draw.size(); i++) draw[i] = draw[i] / total;
  return draw;
}

// rbinom(1, size, prob): a whole number, or a missing value.
inline double rbinom_one(double size, double prob) {
  double draw = R::rbinom(size, prob);
  if (ISNAN(draw) || draw > INT_MAX || draw <= INT_MIN) return NA_REAL;
  return (double) (int) draw;
}

// One draw from a multinomial distribution, as a chain of binomial draws.
inline dvec rmultinomial(double size, const dvec& prob) {
  std::size_t K = prob.size();
  dvec draw(K, 0.0);
  double remaining = r_sum(prob);
  for (std::size_t i = 0; i < K; i++) {
    if (ISNAN(prob[i])) stop_missing_condition();
    if (prob[i] != 0) {
      double p = prob[i] / remaining;
      if (ISNAN(p)) stop_missing_condition();
      draw[i] = p < 1 ? rbinom_one(size, p) : size;
      size = size - draw[i];
    }
    if (ISNAN(size)) stop_missing_condition();
    if (size <= 0) return draw;
    remaining = remaining - prob[i];
  }
  if (K >= 2) draw[K - 2] = size;  // as the R code: the count left over goes to category K - 1
  return draw;
}

// Parametric alternative f1 shifted to start right above the cutoff C.
inline dvec shifted_density(Model model, const dvec& x, double xi, double delta, double nu,
                            double C) {
  dvec f1(x.size(), 0.0);
  for (std::size_t i = 0; i < x.size(); i++) {
    if (!(x[i] > C)) continue;
    double shift = x[i] - C - 1;
    if (model == MODEL_P) {
      f1[i] = r_exp(shift * r_log(xi) - r_lgamma(shift + 1) - xi);
    } else if (model == MODEL_GP) {
      f1[i] = r_exp(r_log(delta) - r_lgamma(shift + 1) +
                    (shift - 1) * r_log(delta + nu * shift) - (delta + nu * shift));
    }
  }
  return f1;
}

}  // namespace localfdr

#endif
