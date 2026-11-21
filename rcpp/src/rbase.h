// -----------------------------------------------------------------------------
// File:    rcpp/src/rbase.h
// Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
//          false discovery rate for sparse count data with application to the
//          discovery of hotspots in protein domains. Ann. Appl. Stat.
//          16(3):1459-1475. Tables 1 to 3, Figures S1 to S3.
// Purpose: What base R supplies to the code in R/, in C++ and with R's order
//          of operations: sum() and prod() added in long double, %% of
//          doubles, log(), exp() and lgamma() that keep a missing value,
//          min(), max(), pmin() and pmax() with R's treatment of missing
//          values, the matrix product, chol() and eigen() of a symmetric
//          matrix through the BLAS and LAPACK library R is linked to, called
//          as R's own C code calls them.
// Inputs:  None.
// Outputs: None.
// Usage:   #include "rbase.h"   (first, in every file of rcpp/src/)
//          Compiled by rcpp/R/load.R with -ffp-contract=off, so that a
//          product and a sum are never fused into one rounding: R rounds
//          after each operation.
// Source:  R 4.6.1: src/main/summary.c (rsum, rprod, min, max, pmin, pmax),
//          src/main/arithmetic.c (myfmod, R_log, math1), src/main/array.c
//          (matprod, simple_matprod), src/modules/lapack/Lapack.c (La_chol,
//          La_rs), base::eigen, base::seq.default.
// -----------------------------------------------------------------------------
#ifndef LOCALFDR_RBASE_H
#define LOCALFDR_RBASE_H

#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
#define FCONE
#endif

#include <cfloat>
#include <cmath>
#include <string>
#include <vector>

#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#endif

namespace localfdr {

typedef std::vector<double> dvec;

// ----- Elements and missing values -----

// x[i] as R reads it: a missing value beyond the end of the vector.
inline double at(const dvec& x, std::size_t i) {
  return i < x.size() ? x[i] : NA_REAL;
}

// anyNA(x).
inline bool any_na(const dvec& x) {
  for (std::size_t i = 0; i < x.size(); i++) {
    if (ISNAN(x[i])) return true;
  }
  return false;
}

// The error of `if (NA)`.
inline void stop_missing_condition() {
  Rcpp::stop("missing value where TRUE/FALSE needed");
}

// ----- Functions of one number: a missing value stays as it is (math1) -----

inline double r_log(double x) {
  if (ISNAN(x)) return x;
  return x > 0 ? std::log(x) : x == 0 ? R_NegInf : R_NaN;
}

inline double r_exp(double x) {
  return ISNAN(x) ? x : std::exp(x);
}

inline double r_lgamma(double x) {
  return ISNAN(x) ? x : R::lgammafn(x);
}

// ----- Sums and products, added in long double as rsum() and rprod() -----

struct Sum {
  long double s;
  Sum() : s(0.0) {}
  void add(double x) { s += x; }
  double value() const {
    if (s > DBL_MAX) return R_PosInf;
    if (s < -DBL_MAX) return R_NegInf;
    double out = 0.0;
    out += (double) s;
    return out;
  }
};

inline double r_sum(const dvec& x) {
  Sum total;
  for (std::size_t i = 0; i < x.size(); i++) total.add(x[i]);
  return total.value();
}

inline double r_prod(const dvec& x) {
  long double s = 1.0;
  for (std::size_t i = 0; i < x.size(); i++) s *= x[i];
  double out = 1.0;
  out *= (double) s;
  return out;
}

// ----- Smallest and largest -----

// min(a, b) of two numbers.
inline double r_min2(double a, double b) {
  if (ISNAN(a)) return (ISNAN(b) && ISNA(b)) ? b : a;
  if (ISNAN(b)) return b;
  return b < a ? b : a;
}

// max(a, b) of two numbers.
inline double r_max2(double a, double b) {
  if (ISNAN(a)) return (ISNAN(b) && ISNA(b)) ? b : a;
  if (ISNAN(b)) return b;
  return b > a ? b : a;
}

// max(x) of a vector without missing values.
inline double r_max(const dvec& x) {
  double s = R_NegInf;
  bool updated = false;
  for (std::size_t i = 0; i < x.size(); i++) {
    if (!updated || x[i] > s) {
      s = x[i];
      updated = true;
    }
  }
  return s;
}

// pmin(bound, x) and pmax(bound, x) at one element; bound is not missing.
inline double r_pmin(double bound, double x) {
  if (ISNAN(x)) return x;
  return x < bound ? x : bound;
}

inline double r_pmax(double bound, double x) {
  if (ISNAN(x)) return x;
  return x > bound ? x : bound;
}

// ----- Remainder and sequences -----

// x1 %% x2 of doubles (myfmod).
inline double r_mod(double x1, double x2) {
  if (x2 == 0.0) return R_NaN;
  const double q_1_eps = 1 / DBL_EPSILON;
  if (std::fabs(x2) > q_1_eps && R_FINITE(x1) && std::fabs(x1) <= std::fabs(x2)) {
    return (std::fabs(x1) == std::fabs(x2)) ? 0
      : ((x1 < 0 && x2 > 0) || (x2 < 0 && x1 > 0)) ? x1 + x2 : x1;
  }
  double q = x1 / x2;
  long double tmp = (long double) x1 - std::floor(q) * (long double) x2;
  return (double) (tmp - floorl(tmp / x2) * x2);
}

// length(seq(1, to, 1)), with the error of seq() when `to` is below 1.
inline double seq_length_from_1(double to) {
  double del = to - 1;
  if (ISNAN(del) || !R_FINITE(del)) Rcpp::stop("'to' must be a finite number");
  if (del < 0) Rcpp::stop("wrong sign in 'by' argument");
  return std::floor(del + 1e-10) + 1;
}

// The whole numbers from:to as R lists them, upwards or downwards.
inline std::vector<double> r_colon(double from, double to) {
  if (ISNAN(from) || ISNAN(to)) Rcpp::stop("NA/NaN argument");
  std::vector<double> out;
  std::size_t n = (std::size_t) (std::fabs(to - from) + 1 + FLT_EPSILON);
  for (std::size_t k = 0; k < n; k++) out.push_back(from <= to ? from + k : from - k);
  return out;
}

// times of rep(x, times): cut to a whole number.
inline std::size_t rep_times(double times) {
  if (ISNAN(times) || times < 0) Rcpp::stop("invalid 'times' argument");
  return (std::size_t) times;
}

// ----- Matrices (column-major, as R stores them) -----

struct Mat {
  int nr, nc;
  dvec v;
  Mat() : nr(0), nc(0) {}
  Mat(int rows, int cols) : nr(rows), nc(cols), v((std::size_t) rows * cols, 0.0) {}
  double& operator()(int i, int j) { return v[i + (std::size_t) nr * j]; }
  double operator()(int i, int j) const { return v[i + (std::size_t) nr * j]; }
};

inline Mat as_mat(SEXP x) {
  Rcpp::NumericVector values(x);
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  Mat out;
  if (Rf_isNull(dim)) {  // as.matrix(): a vector is one column
    out.nr = (int) values.size();
    out.nc = 1;
  } else {
    if (Rf_length(dim) != 2) Rcpp::stop("a matrix is needed");
    out.nr = INTEGER(dim)[0];
    out.nc = INTEGER(dim)[1];
  }
  out.v.assign(values.begin(), values.end());
  return out;
}

inline Rcpp::NumericMatrix as_r(const Mat& m) {
  Rcpp::NumericMatrix out(m.nr, m.nc);
  std::copy(m.v.begin(), m.v.end(), out.begin());
  return out;
}

inline Mat transpose(const Mat& m) {
  Mat out(m.nc, m.nr);
  for (int i = 0; i < m.nr; i++) {
    for (int j = 0; j < m.nc; j++) out(j, i) = m(i, j);
  }
  return out;
}

// The test of matprod() for values that are not finite; it adds neighbours,
// as R does.
inline bool may_have_nan_or_inf(const dvec& x) {
  std::size_t n = x.size();
  if ((n & 1) != 0 && !R_FINITE(x[0])) return true;
  for (std::size_t i = n & 1; i < n; i += 2) {
    if (!R_FINITE(x[i] + x[i + 1])) return true;
  }
  return false;
}

// x %*% y: dgemv when y has one column or x one row, dgemm otherwise; summed
// in a loop when a value is not finite.
inline Mat matprod(const Mat& x, const Mat& y) {
  if (x.nc != y.nr) Rcpp::stop("non-conformable arguments");
  Mat z(x.nr, y.nc);
  int nrx = x.nr, ncx = x.nc, nry = y.nr, ncy = y.nc;
  if (nrx == 0 || ncx == 0 || nry == 0 || ncy == 0) return z;
  if (may_have_nan_or_inf(x.v) || may_have_nan_or_inf(y.v)) {
    for (int i = 0; i < nrx; i++) {
      for (int k = 0; k < ncy; k++) {
        long double sum = 0.0;
        for (int j = 0; j < ncx; j++) sum += x(i, j) * y(j, k);
        z(i, k) = (double) sum;
      }
    }
    return z;
  }
  double one = 1.0, zero = 0.0;
  int ione = 1;
  if (ncy == 1) {
    F77_CALL(dgemv)("N", &nrx, &ncx, &one, x.v.data(), &nrx, y.v.data(), &ione, &zero,
                    z.v.data(), &ione FCONE);
  } else if (nrx == 1) {
    F77_CALL(dgemv)("T", &nry, &ncy, &one, y.v.data(), &nry, x.v.data(), &ione, &zero,
                    z.v.data(), &ione FCONE);
  } else {
    F77_CALL(dgemm)("N", "N", &nrx, &ncy, &ncx, &one, x.v.data(), &nrx, y.v.data(), &nry,
                    &zero, z.v.data(), &nrx FCONE FCONE);
  }
  return z;
}

// chol(A): the upper triangular factor, by LAPACK dpotrf on the upper
// triangle (La_chol).
inline Mat chol(const Mat& A) {
  if (A.nr != A.nc) Rcpp::stop("'a' must be a square matrix");
  if (A.nr <= 0) Rcpp::stop("'a' must have dims > 0");
  Mat ans = A;
  int n = A.nr, info = 0;
  for (int j = 0; j < n; j++) {
    for (int i = j + 1; i < n; i++) ans(i, j) = 0.0;
  }
  F77_CALL(dpotrf)("U", &n, ans.v.data(), &n, &info FCONE);
  if (info > 0) Rcpp::stop("the leading minor of order %d is not positive", info);
  if (info < 0) Rcpp::stop("argument %d of Lapack routine %s had invalid value", -info, "dpotrf");
  return ans;
}

inline bool exactly_symmetric(const Mat& M) {
  for (int i = 0; i < M.nr; i++) {
    for (int j = 0; j < i; j++) {
      if (!(M(i, j) == M(j, i))) return false;
    }
  }
  return true;
}

// eigen(M): values, largest first, and vectors in the columns of `vectors`.
// A symmetric matrix goes to LAPACK dsyevr as R calls it (La_rs). A matrix
// that is not exactly symmetric is passed to R's eigen(), which decides with
// a tolerance whether to treat it as symmetric.
inline void eigen(const Mat& M, dvec& values, Mat& vectors) {
  int n = M.nr;
  if (n == 0) Rcpp::stop("0 x 0 matrix");
  if (n != M.nc) Rcpp::stop("non-square matrix in 'eigen'");
  for (std::size_t i = 0; i < M.v.size(); i++) {
    if (!R_FINITE(M.v[i])) Rcpp::stop("infinite or missing values in 'x'");
  }
  values.assign(n, 0.0);
  vectors = Mat(n, n);
  if (!exactly_symmetric(M)) {
    Rcpp::Function r_eigen("eigen", R_BaseEnv);
    Rcpp::List e = r_eigen(as_r(M));
    if (TYPEOF(e["values"]) != REALSXP || TYPEOF(e["vectors"]) != REALSXP) {
      Rcpp::stop("the eigenvalues of the matrix are not real");
    }
    Rcpp::NumericVector r_values = e["values"];
    Rcpp::NumericVector r_vectors = e["vectors"];
    std::copy(r_values.begin(), r_values.end(), values.begin());
    std::copy(r_vectors.begin(), r_vectors.end(), vectors.v.begin());
    return;
  }
  Mat a = M, z(n, n);
  dvec w(n);
  std::vector<int> isuppz(2 * (std::size_t) n);
  double vl = 0.0, vu = 0.0, abstol = 0.0, work_size = 0.0;
  int il = 0, iu = 0, m = 0, info = 0, lwork = -1, liwork = -1, iwork_size = 0;
  F77_CALL(dsyevr)("V", "A", "L", &n, a.v.data(), &n, &vl, &vu, &il, &iu, &abstol, &m, w.data(),
                   z.v.data(), &n, isuppz.data(), &work_size, &lwork, &iwork_size, &liwork,
                   &info FCONE FCONE FCONE);
  if (info != 0) Rcpp::stop("error code %d from Lapack routine '%s'", info, "dsyevr");
  lwork = (int) work_size;
  liwork = iwork_size;
  dvec work(lwork);
  std::vector<int> iwork(liwork);
  F77_CALL(dsyevr)("V", "A", "L", &n, a.v.data(), &n, &vl, &vu, &il, &iu, &abstol, &m, w.data(),
                   z.v.data(), &n, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork,
                   &info FCONE FCONE FCONE);
  if (info != 0) Rcpp::stop("error code %d from Lapack routine '%s'", info, "dsyevr");
  for (int j = 0; j < n; j++) {
    values[j] = w[n - 1 - j];
    for (int i = 0; i < n; i++) vectors(i, j) = z(i, n - 1 - j);
  }
}

}  // namespace localfdr

#endif
