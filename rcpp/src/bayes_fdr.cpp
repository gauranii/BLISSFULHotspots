// -----------------------------------------------------------------------------
// File:    rcpp/src/bayes_fdr.cpp
// Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
//          false discovery rate for sparse count data with application to the
//          discovery of hotspots in protein domains. Ann. Appl. Stat.
//          16(3):1459-1475. Tables 1 to 3, Figures S1 and S3.
// Purpose: The Gibbs sampler of R/bayes_fdr.R in C++: bayes_local_fdr() and
//          every function it calls, as R functions of the same names and
//          arguments. The whole chain runs in C++, for the four alternatives
//          (NP, SP, P, GP), the three versions of the sampler
//          (SAMPLER_VERSIONS, read from R) and the zigp and cutoff switches.
//          Every step follows the R code operation by operation and draws
//          from R's random number generator in the same order, so a seeded
//          chain is the chain R gives.
// Inputs:  None.
// Outputs: None.
// Usage:   source("R/setup.R"); source("rcpp/R/load.R")
//          fit <- bayes_local_fdr(x, model = "NP", gamma = 1)
//          Two steps call R: the covariance of the chain behind the adapted
//          proposal (stats::var, once every 100 iterations; see
//          adapted_proposal() below) and the two functions stored in a
//          sampler version (extra_atoms, fdr_from).
// Source:  R/bayes_fdr.R; R 4.6.1 src/library/stats/src/random.c (how the
//          random number functions call the generators), src/main/array.c
//          (colMeans).
// -----------------------------------------------------------------------------
// [[Rcpp::plugins(localfdr)]]
#include "distributions_core.h"

using namespace localfdr;

namespace {

// ----- What the R code passes around in lists -----

// An entry of SAMPLER_VERSIONS.
struct Version {
  bool check_range, replace_missing, no_cutoff_one, sp_flat;
  double logit_scale, theta_max;
  SEXP extra_atoms, fdr_from;
};

template <typename T>
T element(const Rcpp::List& list, const char* name, T otherwise) {
  if (!list.containsElementNamed(name) || Rf_isNull(list[name])) return otherwise;
  return Rcpp::as<T>(list[name]);
}

dvec vector_element(const Rcpp::List& list, const char* name) {
  if (!list.containsElementNamed(name) || Rf_isNull(list[name])) return dvec();
  return Rcpp::as<dvec>(list[name]);
}

Version read_version(const Rcpp::List& version) {
  Version v;
  v.check_range = element<bool>(version, "check_range", false);
  v.replace_missing = element<bool>(version, "replace_missing", false);
  v.no_cutoff_one = element<std::string>(version, "no_cutoff", "") == "one";
  v.sp_flat = element<std::string>(version, "sp_prior", "") == "flat";
  v.logit_scale = element<double>(version, "logit_scale", NA_REAL);
  v.theta_max = element<double>(version, "theta_max", NA_REAL);
  v.extra_atoms = version.containsElementNamed("extra_atoms") ? (SEXP) version["extra_atoms"]
                                                              : R_NilValue;
  v.fdr_from = version.containsElementNamed("fdr_from") ? (SEXP) version["fdr_from"] : R_NilValue;
  return v;
}

// State of the chain.
struct State {
  double eta, lambda, theta, pi0, xi, delta, nu, tau;
  dvec psi_sp, psi_np;
};

State read_state(const Rcpp::List& s) {
  State out;
  out.eta = element<double>(s, "eta", NA_REAL);
  out.lambda = element<double>(s, "lambda", NA_REAL);
  out.theta = element<double>(s, "theta", NA_REAL);
  out.pi0 = element<double>(s, "pi0", NA_REAL);
  out.xi = element<double>(s, "xi", NA_REAL);
  out.delta = element<double>(s, "delta", NA_REAL);
  out.nu = element<double>(s, "nu", NA_REAL);
  out.tau = element<double>(s, "tau", NA_REAL);
  out.psi_sp = vector_element(s, "psi_sp");
  out.psi_np = vector_element(s, "psi_np");
  return out;
}

// The cutoff C, or NULL when none could be drawn.
struct Cutoff {
  bool null;
  double value;
  std::size_t index;  // place among the distinct counts, when drawn
};

Cutoff read_cutoff(SEXP C) {
  Cutoff out = {true, NA_REAL, 0};
  if (!Rf_isNull(C)) {
    out.null = false;
    out.value = Rcpp::as<double>(C);
  }
  return out;
}

Cutoff cutoff_at(double value) {
  Cutoff out = {false, value, 0};
  return out;
}

// Output of draw_groups().
struct Groups {
  double K, N0, N1, N2, N3, n00;
  dvec n0, n1;
};

Groups read_groups(const Rcpp::List& groups) {
  Groups g;
  g.K = element<double>(groups, "K", NA_REAL);
  g.N0 = element<double>(groups, "N0", NA_REAL);
  g.N1 = element<double>(groups, "N1", NA_REAL);
  g.N2 = element<double>(groups, "N2", NA_REAL);
  g.N3 = element<double>(groups, "N3", NA_REAL);
  g.n00 = element<double>(groups, "n00", NA_REAL);
  g.n0 = vector_element(groups, "n0");
  g.n1 = vector_element(groups, "n1");
  return g;
}

// A Gamma prior (shape kappa, rate zeta).
struct GammaPrior {
  double kappa, zeta;
};

GammaPrior read_gamma_prior(const Rcpp::List& prior) {
  GammaPrior out = {element<double>(prior, "kappa", NA_REAL), element<double>(prior, "zeta", NA_REAL)};
  return out;
}

// One element of a vector, of the type of the vector.
SEXP element_of(SEXP x, std::size_t i) {
  if (TYPEOF(x) == INTSXP) return Rcpp::wrap(INTEGER(x)[i]);
  double value = Rcpp::as<dvec>(x)[i];
  return Rcpp::wrap(value);
}

// ----- Mixture density and local fdr -----

// pi0 f0(u): the default of the argument f0 of mixture_fdr().
dvec null_density(const dvec& u, const State& s, const Version& v) {
  dvec f0 = dzigp(u, s.eta, s.lambda, s.theta, v.check_range);
  for (std::size_t i = 0; i < f0.size(); i++) f0[i] = s.pi0 * f0[i];
  return f0;
}

// The mixture density f at the distinct counts.
void mixture_density(Model model, const dvec& u, const State& s, const Cutoff& C,
                     const Version& v, const dvec& f0, dvec& f) {
  double cutoff = C.value;
  if (C.null) {
    if (!v.no_cutoff_one) Rcpp::stop("no cutoff could be drawn");
    cutoff = 1;
  }
  dvec f1(u.size(), 0.0);
  if (model == MODEL_P || model == MODEL_GP) {
    f1 = shifted_density(model, u, s.xi, s.delta, s.nu, cutoff);
    for (std::size_t i = 0; i < f1.size(); i++) f1[i] = (1 - s.pi0) * f1[i];
  } else if (model == MODEL_SP) {
    f1 = s.psi_sp;
    for (std::size_t i = 0; i < f1.size(); i++) f1[i] = (1 - s.pi0) * f1[i];
  }
  f.assign(u.size(), 0.0);
  for (std::size_t i = 0; i < u.size(); i++) {
    if (u[i] <= cutoff) {
      f[i] = at(f0, i);
    } else if (u[i] > cutoff) {
      f[i] = model != MODEL_NP ? at(f0, i) + at(f1, i) : at(s.psi_np, i);
    }
  }
}

// The local fdr min(1, pi0 f0 / f).
dvec local_fdr(const dvec& f0, const dvec& f) {
  dvec fdr(f.size());
  for (std::size_t i = 0; i < f.size(); i++) fdr[i] = r_pmin(1, at(f0, i) / f[i]);
  return fdr;
}

dvec mixture_fdr_values(Model model, const dvec& u, const State& s, const Cutoff& C,
                        const Version& v) {
  if (C.null && !v.no_cutoff_one) Rcpp::stop("no cutoff could be drawn");
  dvec f0 = null_density(u, s, v), f;
  mixture_density(model, u, s, C, v, f0, f);
  return local_fdr(f0, f);
}

// ----- Full conditional draws -----

Cutoff draw_cutoff_c(Model model, const dvec& u, const dvec& frequency, const State& s,
                     const GammaPrior& tau_prior, const Version& v) {
  const Cutoff none = {true, NA_REAL, 0};
  double kappa = tau_prior.kappa, zeta = tau_prior.zeta;
  std::size_t U = u.size();
  if (any_na(s.psi_sp) || any_na(s.psi_np)) Rcpp::stop("an atom probability is missing");
  dvec loglik(U, 0.0), log_poisson(U), log_f(U), f;
  for (std::size_t k = 0; k < U; k++) {
    log_poisson[k] = u[k] * r_log(s.tau) - r_lgamma(u[k] + 1) - s.tau;
  }
  double log_gamma = kappa * r_log(zeta) - r_lgamma(kappa) + (kappa - 1) * r_log(s.tau) -
    (s.tau * zeta);
  dvec f0 = null_density(u, s, v);

  for (std::size_t k = 0; k < U; k++) {
    mixture_density(model, u, s, cutoff_at(u[k]), v, f0, f);
    bool not_a_number = false;
    for (std::size_t i = 0; i < U; i++) {
      log_f[i] = r_log(f[i]);
      if (R_IsNaN(log_f[i])) not_a_number = true;
    }
    if (!not_a_number) {
      Sum total;
      for (std::size_t i = 0; i < U; i++) total.add(at(frequency, i) * log_f[i]);
      loglik[k] = total.value();
    }
  }

  dvec logpost(U);
  for (std::size_t k = 0; k < U; k++) {
    logpost[k] = loglik[k] + log_poisson[k] + log_gamma;
    if (ISNAN(logpost[k])) return none;
  }
  double largest = r_max(logpost);
  dvec weight(U);
  for (std::size_t k = 0; k < U; k++) weight[k] = r_exp(logpost[k] - largest);
  double total = r_sum(weight);
  dvec q(U);
  for (std::size_t k = 0; k < U; k++) {
    q[k] = weight[k] / total;
    if (ISNAN(q[k])) return none;
  }
  dvec draw = rmultinomial(1, q);
  std::size_t which = 0;  // which.max(): the first of the largest
  for (std::size_t k = 1; k < U; k++) {
    if (draw[k] > draw[which]) which = k;
  }
  Cutoff out = {false, u[which], which};
  return out;
}

double draw_tau_c(const Cutoff& C, const GammaPrior& tau_prior) {
  double shape = C.null ? tau_prior.kappa : C.value + tau_prior.kappa;
  double rate = tau_prior.zeta + 1;
  return R::rgamma(shape, 1 / rate);
}

Groups draw_groups_c(Model model, const dvec& u, const dvec& x_sorted, const dvec& frequency,
                     const State& s, const Cutoff& C, const Version& v) {
  std::size_t N = x_sorted.size();
  if (N == 0) Rcpp::stop("no observations");
  double K = x_sorted[N - 1];
  if (ISNAN(K) || K < 0) Rcpp::stop("the counts must be non-negative");
  std::size_t bins = (std::size_t) K + 1;  // counts 0, ..., K
  double cutoff = C.null ? (v.no_cutoff_one ? 1 : K) : C.value;

  dvec fdr = mixture_fdr_values(model, u, s, cutoff_at(cutoff), v);
  dvec p;
  p.reserve(N);
  for (std::size_t i = 0; i < fdr.size(); i++) {
    std::size_t times = rep_times(at(frequency, i));
    for (std::size_t r = 0; r < times; r++) p.push_back(fdr[i]);
  }
  if (p.size() != N) Rcpp::stop("the frequencies do not add up to the number of observations");
  for (std::size_t i = 0; i < N; i++) p[i] = r_pmax(x_sorted[i] <= cutoff ? 1 : 0, p[i]);
  if (v.replace_missing) {
    for (std::size_t i = 0; i < N; i++) {
      if (ISNAN(p[i])) p[i] = R::runif(0.0, 1.0);
    }
  }
  dvec z(N);
  bool missing = false;
  for (std::size_t i = 0; i < N; i++) {
    z[i] = rbinom_one(1, p[i]);
    if (ISNAN(z[i])) missing = true;
  }

  dvec n0(bins, 0.0), n1(bins, 0.0), all(bins, 0.0);
  if (missing) {
    for (std::size_t value = 0; value < bins; value++) {
      Sum null, alternative;
      for (std::size_t i = 0; i < N; i++) {
        double here = x_sorted[i] == (double) value ? 1 : 0;
        null.add(z[i] * here);
        alternative.add((1 - z[i]) * here);
      }
      n0[value] = null.value();
      n1[value] = alternative.value();
    }
  }
  for (std::size_t i = 0; i < N; i++) {  // tabulate(x + 1, K + 1)
    double bin = std::trunc(x_sorted[i] + 1);
    if (!(bin >= 1 && bin <= (double) bins)) continue;
    std::size_t place = (std::size_t) bin - 1;
    all[place] += 1;
    if (!missing) {
      if (z[i] == 1) n0[place] += 1;
      if (z[i] == 0) n1[place] += 1;
    }
  }

  Groups g;
  Sum N0, N1, N2, N3;
  for (std::size_t value = 0; value < bins; value++) {
    double above = (double) value > cutoff ? 1 : 0;
    N0.add(n0[value]);
    N1.add(n1[value] * above);
    N2.add((double) value * n0[value]);
    N3.add((double) value * n1[value] * above);
  }
  g.K = K;
  g.N0 = N0.value();
  g.N1 = N1.value();
  g.N2 = N2.value();
  g.N3 = N3.value();
  g.n00 = n0[0];
  for (std::size_t value = 0; value < bins; value++) {
    if (all[value] > 0) {
      g.n0.push_back(n0[value]);
      g.n1.push_back(n1[value]);
    }
  }
  return g;
}

double draw_xi_c(const Cutoff& C, const Groups& groups, const GammaPrior& xi_prior) {
  double cutoff = C.null ? groups.K : C.value;
  double shape = groups.N3 + xi_prior.kappa - (cutoff + 1) * groups.N1 + 1;
  double rate = groups.N1 + xi_prior.zeta;
  return R::rgamma(shape, 1 / rate);
}

dvec concentration_above_cutoff_c(double C, double P, double gamma, double epsilon) {
  if (ISNAN(C)) stop_missing_condition();
  std::size_t atoms = rep_times(P);
  if (atoms == 0) Rcpp::stop("the number of atoms must be positive");
  if (C == -1) return dvec(atoms, gamma);
  double n = seq_length_from_1(C + 1);
  dvec out(atoms);
  for (std::size_t i = 0; i < atoms; i++) out[i] = (double) i < n ? epsilon : gamma;
  return out;
}

dvec draw_atoms_c(const dvec& u, const dvec& psi, const dvec& concentration, const State& s,
                  double gamma, const Cutoff& C, double atoms, const Version& v) {
  double cutoff = C.value;
  if (C.null) {
    if (!v.no_cutoff_one) Rcpp::stop("no cutoff could be drawn");
    cutoff = 0;
  }
  std::size_t U = u.size();
  dvec u_below, alternative;
  Sum below;
  for (std::size_t i = 0; i < U; i++) {
    if (u[i] <= cutoff) {
      u_below.push_back(u[i]);
      below.add(at(psi, i));
    } else if (u[i] > cutoff) {
      alternative.push_back(at(concentration, i));
    }
  }
  double null_mass = r_min2(1 - 1e-10, below.value());
  if (atoms < (double) U) {
    if (Rf_isNull(v.extra_atoms)) Rcpp::stop("the sampler version has no extra_atoms");
    Rcpp::Function extra_atoms(v.extra_atoms);
    atoms = Rcpp::as<double>(extra_atoms((int) U, r_max(u)));
  }
  std::size_t extra = rep_times(r_max2(0, atoms - (double) U + cutoff + 1));
  alternative.insert(alternative.end(), extra, gamma);
  dvec psi_null = dzigp(u_below, s.eta, s.lambda, s.theta, v.check_range);
  for (std::size_t i = 0; i < psi_null.size(); i++) psi_null[i] = s.pi0 * psi_null[i];
  dvec psi_alternative = rdirichlet(alternative, v.replace_missing);
  for (std::size_t i = 0; i < psi_alternative.size(); i++) {
    psi_alternative[i] = (1 - null_mass) * psi_alternative[i];
  }
  dvec out(U, NA_REAL);
  for (std::size_t i = 0; i < U; i++) {
    if (i < psi_null.size()) {
      out[i] = psi_null[i];
    } else if (i - psi_null.size() < psi_alternative.size()) {
      out[i] = psi_alternative[i - psi_null.size()];
    }
  }
  return out;
}

// ----- Random walk Metropolis-Hastings steps -----

dvec null_parameters_c(const dvec& pre, const Version& v) {
  dvec par(3);
  par[0] = R::plogis(at(pre, 0), 0.0, 1.0, 1, 0);
  par[1] = r_exp(at(pre, 1));
  par[2] = R::plogis(at(pre, 2), 0.0, v.logit_scale, 1, 0);
  return par;
}

dvec gp_parameters_c(const dvec& pre, const Version& v) {
  dvec par(2);
  par[0] = r_exp(at(pre, 0));
  par[1] = R::plogis(at(pre, 1), 0.0, v.logit_scale, 1, 0);
  return par;
}

// The two sums over the distinct counts shared by the log likelihoods:
// sum(n[2:last] * (y - 1) * log(rate + dispersion * y)) and
// sum(n[2:last] * lgamma(y + 1)), y = 1, ..., last - 1.
void count_sums(const dvec& n, double rate, double dispersion, double& powers, double& factorials) {
  double length = seq_length_from_1((double) n.size() - 1);
  Sum first, second;
  for (double y = 1; y <= length; y++) {
    double count = at(n, (std::size_t) y);
    first.add(count * (y - 1) * r_log(rate + dispersion * y));
    second.add(count * r_lgamma(y + 1));
  }
  powers = first.value();
  factorials = second.value();
}

double loglik_null_c(const dvec& par, const Groups& groups) {
  double eta = at(par, 0), lambda = at(par, 1), theta = at(par, 2);
  double powers, factorials;
  count_sums(groups.n0, lambda, theta, powers, factorials);
  return groups.n00 * r_log(eta + (1 - eta) * r_exp(-lambda)) +
    (groups.N0 - groups.n00) * r_log(1 - eta) +
    (groups.N0 - groups.n00) * (r_log(lambda) - lambda) - groups.N2 * theta +
    powers - factorials;
}

double loglik_gp_c(const dvec& par, const Groups& groups) {
  double delta = at(par, 0), nu = at(par, 1);
  double powers, factorials;
  count_sums(groups.n1, delta, nu, powers, factorials);
  return groups.N1 * r_log(delta) - delta * groups.N1 - nu * groups.N3 + powers - factorials;
}

double logpost_null_c(const dvec& pre, const Groups& groups, const Version& v) {
  dvec par = null_parameters_c(pre, v);
  double prior = R::dunif(par[0], 0.0, 1.0, 1) + -0.5 * r_log(par[1]) +
    R::dunif(par[2], 0.0, v.theta_max, 1);
  dvec jacobian(3);
  jacobian[0] = R::dlogis(at(pre, 0), 0.0, 1.0, 0);
  jacobian[1] = r_exp(at(pre, 1));
  jacobian[2] = R::dlogis(at(pre, 2), 0.0, v.logit_scale, 0);
  return loglik_null_c(par, groups) + prior + r_log(std::fabs(r_prod(jacobian)));
}

double logpost_gp_c(const dvec& pre, const Groups& groups, const Version& v) {
  dvec par = gp_parameters_c(pre, v);
  double prior = -0.5 * r_log(par[0]) + R::dunif(par[1], 0.0, 1.0, 1);
  dvec jacobian(2);
  jacobian[0] = r_exp(at(pre, 0));
  jacobian[1] = R::dlogis(at(pre, 1), 0.0, v.logit_scale, 0);
  return loglik_gp_c(par, groups) + prior + r_log(std::fabs(r_prod(jacobian)));
}

Mat positive_definite_c(const Mat& M) {
  int d = M.nr;
  if (M.nc != d) Rcpp::stop("Input matrix is not square!");
  dvec values;
  Mat vectors;
  eigen(M, values, vectors);
  double largest = R_NegInf;
  for (int i = 0; i < d; i++) largest = r_max2(largest, std::fabs(values[i]));
  double tolerance = d * largest * DBL_EPSILON;
  Mat lift(d, d);
  for (int i = 0; i < d; i++) lift(i, i) = r_pmax(0, 2 * tolerance - values[i]);
  Mat correction = matprod(matprod(vectors, lift), transpose(vectors));
  Mat out = M;
  for (std::size_t i = 0; i < out.v.size(); i++) out.v[i] = M.v[i] + correction.v[i];
  return out;
}

// The chain so far, one row per iteration.
struct History {
  int columns;
  dvec values;  // row after row
  std::size_t rows() const { return values.size() / columns; }
};

// Cholesky factor of the adapted proposal covariance.
//
// The covariance of the rows of the chain is computed by R's var(), called
// from here. R adds the products of deviations in a compiled loop whose
// rounding depends on how the compiler built R (it may fuse a product and a
// sum into one rounding), so the same loop compiled here could differ from R
// in the last digit. The step runs once every 100 iterations, and everything
// after it (the eigenvalues, the lift and the Cholesky factor) is done here.
Mat adapted_proposal_c(const History& history, double t, double prop, const Mat& fallback) {
  std::vector<double> rows = r_colon(std::floor(t * prop), t);
  std::vector<std::size_t> kept;
  for (std::size_t i = 0; i < rows.size(); i++) {
    if (rows[i] < 0 || rows[i] > (double) history.rows()) Rcpp::stop("subscript out of bounds");
    if (rows[i] != 0) kept.push_back((std::size_t) rows[i] - 1);
  }
  double length = (double) rows.size();
  Rcpp::Function var("var", Rcpp::Environment::namespace_env("stats"));
  Mat sigma;
  if (kept.size() == 1) {  // one row is a vector in R: its variance is one number
    Rcpp::NumericVector row(history.columns);
    for (int j = 0; j < history.columns; j++) row[j] = history.values[kept[0] * history.columns + j];
    sigma = as_mat(var(row));
  } else {
    Rcpp::NumericMatrix chain((int) kept.size(), history.columns);
    for (std::size_t i = 0; i < kept.size(); i++) {
      for (int j = 0; j < history.columns; j++) {
        chain((int) i, j) = history.values[kept[i] * history.columns + j];
      }
    }
    sigma = as_mat(var(chain));
  }
  for (std::size_t i = 0; i < sigma.v.size(); i++) sigma.v[i] = (length - 1) * sigma.v[i] / length;
  try {
    return chol(positive_definite_c(sigma));
  } catch (Rcpp::exception&) {
    return chol(positive_definite_c(fallback));
  }
}

// One adaptive random walk Metropolis-Hastings step of d parameters.
struct Walk {
  int d;
  dvec pre;       // current point, unconstrained scale
  Mat half;       // Cholesky factor of the proposal covariance
  Mat start_var;  // starting proposal covariance
  History history;
};

// candidate <- t(half) %*% rnorm(d) + pre
dvec propose(const Walk& walk) {
  Mat normal(walk.d, 1);
  for (int i = 0; i < walk.d; i++) normal.v[i] = R::rnorm(0.0, 1.0);
  Mat step = matprod(transpose(walk.half), normal);
  if (step.nr != walk.d) Rcpp::stop("non-conformable arrays");
  dvec candidate(walk.d);
  for (int i = 0; i < walk.d; i++) candidate[i] = step.v[i] + walk.pre[i];
  return candidate;
}

// ----- The sampler -----

// match.arg(model, c("NP", "SP", "P", "GP")).
std::string match_model(SEXP model) {
  const char* choices[] = {"NP", "SP", "P", "GP"};
  if (Rf_isNull(model)) return choices[0];
  if (TYPEOF(model) != STRSXP) Rcpp::stop("'arg' must be NULL or a character vector");
  std::vector<std::string> given = Rcpp::as<std::vector<std::string> >(model);
  if (given.size() == 4 && given[0] == "NP" && given[1] == "SP" && given[2] == "P" &&
      given[3] == "GP") {
    return choices[0];
  }
  if (given.size() != 1) Rcpp::stop("'arg' must be of length 1");
  int found = -1, partial = 0;
  for (int i = 0; i < 4; i++) {
    std::string choice = choices[i];
    if (choice == given[0]) return choice;
    if (!given[0].empty() && choice.compare(0, given[0].size(), given[0]) == 0) {
      found = i;
      partial++;
    }
  }
  if (partial != 1) {
    Rcpp::stop("'arg' should be one of \"NP\", \"SP\", \"P\", \"GP\"");
  }
  return choices[found];
}

}  // namespace

// ----- R functions, with the names and arguments of R/bayes_fdr.R -----

// [[Rcpp::export]]
Rcpp::List mixture_fdr(std::string model, Rcpp::NumericVector u, Rcpp::List s, SEXP C,
                       Rcpp::List version, Rcpp::NumericVector f0) {
  Version v = read_version(version);
  dvec f, null = Rcpp::as<dvec>(f0);
  mixture_density(model_code(model), Rcpp::as<dvec>(u), read_state(s), read_cutoff(C), v, null, f);
  return Rcpp::List::create(Rcpp::Named("f0") = f0, Rcpp::Named("f") = f,
                            Rcpp::Named("fdr") = local_fdr(null, f));
}

// [[Rcpp::export]]
SEXP draw_cutoff(std::string model, SEXP u, Rcpp::NumericVector frequency, Rcpp::List s,
                 Rcpp::List tau_prior, Rcpp::List version) {
  Cutoff C = draw_cutoff_c(model_code(model), Rcpp::as<dvec>(u), Rcpp::as<dvec>(frequency),
                           read_state(s), read_gamma_prior(tau_prior), read_version(version));
  if (C.null) return R_NilValue;
  return element_of(u, C.index);
}

// [[Rcpp::export]]
double draw_tau(SEXP C, Rcpp::List tau_prior) {
  return draw_tau_c(read_cutoff(C), read_gamma_prior(tau_prior));
}

// [[Rcpp::export]]
Rcpp::List draw_groups(std::string model, Rcpp::NumericVector u, SEXP x_sorted,
                       Rcpp::NumericVector frequency, Rcpp::List s, SEXP C, Rcpp::List version) {
  Groups g = draw_groups_c(model_code(model), Rcpp::as<dvec>(u), Rcpp::as<dvec>(x_sorted),
                           Rcpp::as<dvec>(frequency), read_state(s), read_cutoff(C),
                           read_version(version));
  return Rcpp::List::create(
    Rcpp::Named("K") = element_of(x_sorted, Rf_xlength(x_sorted) - 1), Rcpp::Named("N0") = g.N0,
    Rcpp::Named("N1") = g.N1, Rcpp::Named("N2") = g.N2, Rcpp::Named("N3") = g.N3,
    Rcpp::Named("n00") = g.n00, Rcpp::Named("n0") = g.n0, Rcpp::Named("n1") = g.n1);
}

// [[Rcpp::export]]
double draw_xi(SEXP C, Rcpp::List groups, Rcpp::List xi_prior) {
  return draw_xi_c(read_cutoff(C), read_groups(groups), read_gamma_prior(xi_prior));
}

// [[Rcpp::export]]
Rcpp::NumericVector concentration_above_cutoff(double C, double P, double gamma,
                                               double epsilon = 1e-6) {
  return Rcpp::wrap(concentration_above_cutoff_c(C, P, gamma, epsilon));
}

// [[Rcpp::export]]
Rcpp::NumericVector draw_atoms(Rcpp::NumericVector u, Rcpp::NumericVector psi,
                               Rcpp::NumericVector concentration, Rcpp::List s, double gamma,
                               SEXP C, double atoms, Rcpp::List version) {
  return Rcpp::wrap(draw_atoms_c(Rcpp::as<dvec>(u), Rcpp::as<dvec>(psi),
                                 Rcpp::as<dvec>(concentration), read_state(s), gamma,
                                 read_cutoff(C), atoms, read_version(version)));
}

// [[Rcpp::export]]
Rcpp::NumericVector null_parameters(Rcpp::NumericVector pre, Rcpp::List version) {
  return Rcpp::wrap(null_parameters_c(Rcpp::as<dvec>(pre), read_version(version)));
}

// [[Rcpp::export]]
Rcpp::NumericVector gp_parameters(Rcpp::NumericVector pre, Rcpp::List version) {
  return Rcpp::wrap(gp_parameters_c(Rcpp::as<dvec>(pre), read_version(version)));
}

// [[Rcpp::export]]
double loglik_null(Rcpp::NumericVector par, Rcpp::List groups) {
  return loglik_null_c(Rcpp::as<dvec>(par), read_groups(groups));
}

// [[Rcpp::export]]
double loglik_gp(Rcpp::NumericVector par, Rcpp::List groups) {
  return loglik_gp_c(Rcpp::as<dvec>(par), read_groups(groups));
}

// [[Rcpp::export]]
double logpost_null(Rcpp::NumericVector pre, Rcpp::List groups, Rcpp::List version) {
  return logpost_null_c(Rcpp::as<dvec>(pre), read_groups(groups), read_version(version));
}

// [[Rcpp::export]]
double logpost_gp(Rcpp::NumericVector pre, Rcpp::List groups, Rcpp::List version) {
  return logpost_gp_c(Rcpp::as<dvec>(pre), read_groups(groups), read_version(version));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix positive_definite(SEXP M) {
  Rcpp::NumericMatrix out = as_r(positive_definite_c(as_mat(M)));
  SEXP names = Rf_getAttrib(M, R_DimNamesSymbol);
  if (!Rf_isNull(names)) Rf_setAttrib(out, R_DimNamesSymbol, names);
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix adapted_proposal(Rcpp::NumericMatrix history, double t, Rcpp::List adapt,
                                     Rcpp::NumericMatrix fallback) {
  History chain;
  chain.columns = history.ncol();
  chain.values.resize((std::size_t) history.nrow() * history.ncol());
  for (int i = 0; i < history.nrow(); i++) {
    for (int j = 0; j < history.ncol(); j++) {
      chain.values[(std::size_t) i * chain.columns + j] = history(i, j);
    }
  }
  return as_r(adapted_proposal_c(chain, t, element<double>(adapt, "prop", NA_REAL),
                                 as_mat(fallback)));
}

// [[Rcpp::export]]
Rcpp::List bayes_local_fdr(SEXP x, SEXP model, double gamma, Rcpp::List settings, double draws,
                           double burn, double thin, bool zigp, bool cutoff, SEXP version) {
  Model alternative = model_code(match_model(model));
  Rcpp::List versions = Rcpp::Environment::global_env()["SAMPLER_VERSIONS"];
  std::string version_name = Rcpp::as<std::string>(version);
  if (!versions.containsElementNamed(version_name.c_str())) Rcpp::stop("subscript out of bounds");
  Version v = read_version(versions[version_name]);
  Rcpp::List init = settings["init"], adapt = settings["adapt"], proposal = settings["proposal"];
  Rcpp::List pi0_prior = settings["pi0_prior"];
  GammaPrior tau_prior = read_gamma_prior(settings["tau_prior"]);
  GammaPrior xi_prior = read_gamma_prior(settings["xi_prior"]);
  double pi0_alpha = element<double>(pi0_prior, "alpha", NA_REAL);
  double pi0_beta = element<double>(pi0_prior, "beta", NA_REAL);
  double adapt_lower = element<double>(adapt, "lower", NA_REAL);
  double adapt_upper = element<double>(adapt, "upper", NA_REAL);
  double adapt_thin = element<double>(adapt, "thin", NA_REAL);
  double adapt_prop = element<double>(adapt, "prop", NA_REAL);

  // The data: sorted observations, distinct counts and their frequencies.
  double iterations = draws * thin + burn;
  if (ISNAN(iterations) || iterations < 1) Rcpp::stop("the chain needs at least one iteration");
  dvec x_sorted = Rcpp::as<dvec>(x);
  if (x_sorted.empty() || any_na(x_sorted)) Rcpp::stop("the counts must be given, none missing");
  std::sort(x_sorted.begin(), x_sorted.end());
  dvec u;
  for (std::size_t i = 0; i < x_sorted.size(); i++) {
    if (i == 0 || x_sorted[i] != x_sorted[i - 1]) u.push_back(x_sorted[i]);
  }
  std::size_t U = u.size();
  double K = x_sorted.back();
  if (x_sorted.front() < 0) Rcpp::stop("the counts must be non-negative");
  std::vector<int> counted((std::size_t) K + 1, 0);  // c(sum(x == 0), tabulate(x, K))
  for (std::size_t i = 0; i < x_sorted.size(); i++) {
    if (x_sorted[i] == 0) {
      counted[0]++;
    } else {
      double bin = std::trunc(x_sorted[i]);
      if (bin >= 1 && bin <= std::trunc(K)) counted[(std::size_t) bin]++;
    }
  }
  std::vector<int> frequency;
  for (std::size_t i = 0; i < counted.size(); i++) {
    if (counted[i] > 0) frequency.push_back(counted[i]);
  }
  dvec frequency_d(frequency.begin(), frequency.end());

  std::size_t stored = (std::size_t) draws, total = (std::size_t) iterations;
  const char* kept_names[] = {"eta", "lambda", "theta", "pi0", "tau", "C", "xi", "delta", "nu"};
  std::vector<Rcpp::NumericVector> kept;
  for (int i = 0; i < 9; i++) kept.push_back(Rcpp::NumericVector(stored));
  Rcpp::NumericVector acceptance(total);
  dvec fdr_draws(stored * U, NA_REAL);  // column after column

  // State of the chain.
  State s;
  s.eta = element<double>(init, "eta", NA_REAL);
  s.lambda = element<double>(init, "lambda", NA_REAL);
  s.theta = zigp ? element<double>(init, "theta", NA_REAL) : 0;
  s.pi0 = element<double>(init, "pi0", NA_REAL);
  s.xi = element<double>(init, "xi", NA_REAL);
  s.delta = element<double>(init, "delta", NA_REAL);
  s.nu = element<double>(init, "nu", NA_REAL);
  s.tau = element<double>(init, "tau", NA_REAL);
  dvec beta = vector_element(init, "beta");
  s.psi_sp.resize(U);
  for (std::size_t i = 0; i < U; i++) s.psi_sp[i] = at(beta, i);
  s.psi_np = s.psi_sp;
  double atoms = r_min2(1.5 * K, element<double>(init, "atoms", NA_REAL));
  if (any_na(s.psi_sp)) s.psi_sp = rdirichlet(dvec(rep_times(atoms), 0.5), v.replace_missing);
  if (any_na(s.psi_np)) s.psi_np = rdirichlet(dvec(rep_times(atoms), 0.5), v.replace_missing);

  Walk null, gp;
  null.d = 3;
  null.pre = vector_element(proposal, "null_start");
  null.start_var = as_mat(proposal["null_var"]);
  null.half = chol(null.start_var);
  null.history.columns = 3;
  null.history.values.assign(total * 3, NA_REAL);
  gp.d = 2;
  gp.pre = vector_element(proposal, "gp_start");
  gp.start_var = as_mat(proposal["gp_var"]);
  gp.half = chol(gp.start_var);
  gp.history.columns = 2;
  gp.history.values.assign(total * 2, NA_REAL);
  if (null.pre.size() != 3 || gp.pre.size() != 2) Rcpp::stop("the starting points need 3 and 2 values");
  double accepted = 0;
  std::size_t h = 0;
  dvec concentration, fdr;

  for (double t = 1; t <= iterations; t++) {
    if (r_mod(t, 1000) == 0) Rcpp::checkUserInterrupt();
    std::size_t row = (std::size_t) t - 1;
    bool store = (r_mod(t, thin) == 0) && (t > burn);
    bool adapt_now = (t > adapt_lower * iterations) && (r_mod(t, adapt_thin) == 0) &&
      (t < adapt_upper * iterations);

    // Cutoff C and its prior mean tau.
    Cutoff C;
    if (cutoff) {
      C = draw_cutoff_c(alternative, u, frequency_d, s, tau_prior, v);
      s.tau = draw_tau_c(C, tau_prior);
    } else {
      C = cutoff_at(-1);
      s.tau = 0;
    }

    // Group indicators and the null proportion pi0.
    Groups groups = draw_groups_c(alternative, u, x_sorted, frequency_d, s, C, v);
    s.pi0 = R::rbeta(groups.N0 + pi0_alpha, groups.N1 + pi0_beta);

    // Posterior concentration of the atoms.
    if (alternative == MODEL_SP) {
      if (v.sp_flat) {
        concentration.assign(U, gamma);
      } else {
        if (C.null) Rcpp::stop("argument is of length zero");
        concentration = concentration_above_cutoff_c(C.value, (double) U, gamma, 1e-2);
      }
      for (std::size_t i = 0; i < U; i++) concentration[i] = concentration[i] + at(groups.n1, i);
    } else if (alternative == MODEL_NP) {
      concentration.assign(U, gamma);
      for (std::size_t i = 0; i < U; i++) {
        concentration[i] = concentration[i] + at(groups.n0, i) + at(groups.n1, i);
      }
    }

    // ZIGP parameters eta, lambda, theta: adaptive random walk.
    for (int j = 0; j < 3; j++) null.history.values[row * 3 + j] = null.pre[j];
    if (adapt_now) null.half = adapted_proposal_c(null.history, t, adapt_prop, null.start_var);
    dvec candidate = propose(null);
    double log_ratio = logpost_null_c(candidate, groups, v) - logpost_null_c(null.pre, groups, v);
    if (!ISNAN(log_ratio)) {
      if (r_log(R::runif(0.0, 1.0)) < log_ratio) {
        null.pre = candidate;
        accepted = accepted + 1;
      }
    }
    dvec par = null_parameters_c(null.pre, v);
    s.eta = par[0];
    s.lambda = r_min2(par[1], K);
    s.theta = zigp ? par[2] : 0;
    acceptance[row] = accepted * 100 / t;

    // Parameters of a parametric f1.
    if (alternative == MODEL_P) {
      s.xi = draw_xi_c(C, groups, xi_prior);
    } else if (alternative == MODEL_GP) {
      for (int j = 0; j < 2; j++) gp.history.values[row * 2 + j] = gp.pre[j];
      if (adapt_now) gp.half = adapted_proposal_c(gp.history, t, adapt_prop, gp.start_var);
      candidate = propose(gp);
      log_ratio = logpost_gp_c(candidate, groups, v) - logpost_gp_c(gp.pre, groups, v);
      if (!ISNAN(log_ratio)) {
        if (r_log(R::runif(0.0, 1.0)) < log_ratio) gp.pre = candidate;
      }
      par = gp_parameters_c(gp.pre, v);
      s.delta = r_min2(par[0], K);
      s.nu = par[1];
    }

    // Atom probabilities of a semiparametric or nonparametric f1.
    if (alternative == MODEL_SP) {
      s.psi_sp = draw_atoms_c(u, s.psi_sp, concentration, s, gamma, C, atoms, v);
    } else if (alternative == MODEL_NP) {
      s.psi_np = draw_atoms_c(u, s.psi_np, concentration, s, gamma, C, atoms, v);
    }

    fdr = mixture_fdr_values(alternative, u, s, C, v);
    if (store) {
      if (C.null) Rcpp::stop("no cutoff could be drawn");
      if (h >= stored) Rcpp::stop("subscript out of bounds");
      double state[] = {s.eta, s.lambda, s.theta, s.pi0, s.tau, C.value, s.xi, s.delta, s.nu};
      for (int i = 0; i < 9; i++) kept[i][h] = state[i];
      for (std::size_t j = 0; j < U; j++) fdr_draws[h + stored * j] = fdr[j];
      h++;
    }
  }

  // colMeans(fdr_draws[version$fdr_from(burn):draws, ]): each column added in
  // long double and divided by the number of rows.
  if (Rf_isNull(v.fdr_from)) Rcpp::stop("the sampler version has no fdr_from");
  Rcpp::Function fdr_from(v.fdr_from);
  std::vector<double> rows = r_colon(Rcpp::as<double>(fdr_from(burn)), draws);
  std::vector<std::size_t> used;
  for (std::size_t i = 0; i < rows.size(); i++) {
    if (rows[i] < 0 || rows[i] > (double) stored) Rcpp::stop("subscript out of bounds");
    if (rows[i] != 0) used.push_back((std::size_t) rows[i] - 1);
  }
  if (used.size() == 1 || U == 1) Rcpp::stop("'x' must be an array of at least two dimensions");
  Rcpp::NumericVector fdr_mean(U);
  for (std::size_t j = 0; j < U; j++) {
    long double sum = 0.0;
    for (std::size_t i = 0; i < used.size(); i++) sum += fdr_draws[used[i] + stored * j];
    sum /= used.size();
    fdr_mean[j] = (double) sum;
  }

  Rcpp::List out(13);
  Rcpp::CharacterVector names(13);
  names[0] = "values";
  names[1] = "frequency";
  names[2] = "fdr";
  if (TYPEOF(x) == INTSXP) {
    out[0] = Rcpp::IntegerVector(u.begin(), u.end());
  } else {
    out[0] = Rcpp::NumericVector(u.begin(), u.end());
  }
  out[1] = Rcpp::IntegerVector(frequency.begin(), frequency.end());
  out[2] = fdr_mean;
  for (int i = 0; i < 9; i++) {
    out[3 + i] = kept[i];
    names[3 + i] = kept_names[i];
  }
  out[12] = acceptance;
  names[12] = "acceptance";
  out.attr("names") = names;
  return out;
}
