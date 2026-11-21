# -----------------------------------------------------------------------------
# File:    python/localfdr/rbase.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Tables 1 to 3, Figures 1, S1 to S3.
# Purpose: What base R, stats and parallel supply to the code in R/, with R's
#          order of operations: matrix products, chol() and eigen() through
#          the BLAS and LAPACK library R uses, sums, means, variances and
#          covariances that add in R's order, quantile(), min() and max()
#          with R's treatment of missing values, tabulate(), read.csv(),
#          write.csv() with numbers written as R writes them, and mclapply().
# Inputs:  The environment variable LOCALFDR_BLAS names the library (OpenBLAS,
#          set by flake.nix to the one R is linked to). Without it the
#          products and decompositions come from numpy. The library may be
#          built with 32-bit or with 64-bit integers; the calls pass the
#          width it reports.
#          R adds its sums and means in the C type long double, which is a
#          double on arm64 macOS; the sums here are added in doubles, so they
#          round as R's do there.
# Outputs: None.
# Usage:   from localfdr.rbase import matprod, chol, eigen_symmetric, r_var
# Source:  R 4.6.1: src/main/array.c (matprod), src/modules/lapack/Lapack.c
#          (La_chol, La_rs), src/main/summary.c (sum, mean, min, max),
#          src/library/stats/src/cov.c (var), stats::quantile (type 7),
#          src/main/format.c (formatReal, scientific), utils::write.csv,
#          parallel::mclapply.
# -----------------------------------------------------------------------------
import ctypes
import math
import multiprocessing
import os

import numpy as np
import pandas as pd
from numba import njit

from .rmath import muladd

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")  # as flake.nix sets it for R
_LIBRARY = os.environ.get("LOCALFDR_BLAS")
_blas = ctypes.CDLL(_LIBRARY) if _LIBRARY else None


def _integer_type(library):
    """The Fortran INTEGER of the library: 64 bits when OpenBLAS was built
    with INTERFACE64 (its configuration then says USE64BITINT), else 32."""
    if library is None or not hasattr(library, "openblas_get_config"):
        return ctypes.c_int
    library.openblas_get_config.restype = ctypes.c_char_p
    library.openblas_get_config.argtypes = []
    wide = b"USE64BITINT" in library.openblas_get_config()
    return ctypes.c_int64 if wide else ctypes.c_int


_int = _integer_type(_blas)
_int_dtype = np.dtype(_int)
_dbl = ctypes.c_double
_ref = ctypes.byref

# Every argument of the Fortran routines is passed by reference; the lengths
# of the character arguments follow at the end, one size_t each, by value.
_chr, _ptr, _len = ctypes.c_char_p, ctypes.c_void_p, ctypes.c_size_t
if _blas is not None:
    for _name, _characters, _others in [("dgemv_", 1, 10), ("dgemm_", 2, 11),
                                        ("dpotrf_", 1, 4), ("dsyevr_", 3, 18)]:
        _routine = getattr(_blas, _name)
        _routine.argtypes = [_chr] * _characters + [_ptr] * _others + [_len] * _characters
        _routine.restype = None


# ----- Matrix products -----

def as_matrix(x):
    """A float matrix in column-major storage; a vector becomes one column.

    Parameters: x (array). Returns: 2-d Fortran-ordered float64 array.
    """
    x = np.asarray(x, dtype=np.float64)
    if x.ndim == 1:
        x = x.reshape(-1, 1)
    return np.asfortranarray(x)


def _pointer(a):
    return a.ctypes.data_as(ctypes.c_void_p)


def _dgemv(trans, a, x):
    """BLAS dgemv: a %*% x (trans "N") or t(a) %*% x (trans "T")."""
    m, n = a.shape
    y = np.empty(m if trans == b"N" else n)
    _blas.dgemv_(trans, _ref(_int(m)), _ref(_int(n)), _ref(_dbl(1.0)), _pointer(a),
                 _ref(_int(m)), _pointer(x), _ref(_int(1)), _ref(_dbl(0.0)), _pointer(y),
                 _ref(_int(1)), ctypes.c_size_t(1))
    return y


def _dgemm(a, b, rows, cols, inner):
    """BLAS dgemm: a %*% b, a rows x cols matrix."""
    c = np.empty((rows, cols), order="F")
    _blas.dgemm_(b"N", b"N", _ref(_int(rows)), _ref(_int(cols)), _ref(_int(inner)),
                 _ref(_dbl(1.0)), _pointer(a), _ref(_int(a.shape[0])), _pointer(b),
                 _ref(_int(b.shape[0])), _ref(_dbl(0.0)), _pointer(c), _ref(_int(rows)),
                 ctypes.c_size_t(1), ctypes.c_size_t(1))
    return c


@njit(cache=True)
def _simple_matprod(x, y):
    """x %*% y with every entry summed in order (R's simple_matprod)."""
    z = np.empty((x.shape[0], y.shape[1]))
    for i in range(x.shape[0]):
        for k in range(y.shape[1]):
            total = 0.0
            for j in range(x.shape[1]):
                total = muladd(x[i, j], y[j, k], total)
            z[i, k] = total
    return z


def matprod(x, y):
    """`x %*% y` as R computes it.

    dgemv when y has one column or x one row, dgemm otherwise; sums in a
    loop when a value is not finite.

    Parameters: x, y (matrices; a vector is one column).
    Returns: Fortran-ordered matrix.
    """
    x, y = as_matrix(x), as_matrix(y)
    if _blas is None:
        return np.asfortranarray(x @ y)
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        return np.asfortranarray(_simple_matprod(x, y))
    if y.shape[1] == 1:
        return _dgemv(b"N", x, y).reshape(-1, 1, order="F")
    if x.shape[0] == 1:
        return _dgemv(b"T", y, x).reshape(1, -1, order="F")
    return _dgemm(x, y, x.shape[0], y.shape[1], x.shape[1])


# ----- Decompositions -----

def chol(A):
    """`chol(A)`: upper triangular factor R with t(R) %*% R = A, by LAPACK
    dpotrf on the upper triangle, as R calls it.

    Parameters: A (symmetric positive definite matrix).
    Returns: Fortran-ordered upper triangular matrix. Raises ValueError when
    a leading minor is not positive.
    """
    A = as_matrix(A)
    n = A.shape[0]
    if A.shape[1] != n:
        raise ValueError("'a' must be a square matrix")
    if _blas is None:
        if not np.isfinite(A).all():
            raise ValueError("the matrix has values that are not finite")
        try:
            return np.asfortranarray(np.linalg.cholesky(np.triu(A) + np.triu(A, 1).T).T)
        except np.linalg.LinAlgError as error:
            raise ValueError("the leading minor is not positive") from error
    a = np.asfortranarray(np.triu(A))
    info = _int(0)
    _blas.dpotrf_(b"U", _ref(_int(n)), _pointer(a), _ref(_int(n)), _ref(info), ctypes.c_size_t(1))
    if info.value > 0:
        raise ValueError(f"the leading minor of order {info.value} is not positive")
    if info.value < 0:
        raise ValueError(f"argument {-info.value} of Lapack routine dpotrf had invalid value")
    return a


def eigen_symmetric(S):
    """`eigen(S)` of a symmetric matrix: LAPACK dsyevr as R calls it.

    Parameters: S (symmetric matrix).
    Returns: dict with `values` (decreasing) and `vectors` (columns). Raises
    ValueError when a value is not finite or S is not symmetric.
    """
    S = as_matrix(S)
    n = S.shape[0]
    if not np.isfinite(S).all():
        raise ValueError("infinite or missing values in 'x'")
    if not np.allclose(S, S.T, rtol=100 * np.finfo(float).eps, atol=0.0):
        raise ValueError("the matrix is not symmetric")
    if _blas is None:
        values, vectors = np.linalg.eigh(S)
        return {"values": values[::-1].copy(), "vectors": np.asfortranarray(vectors[:, ::-1])}
    a = S.copy(order="F")
    values = np.empty(n)
    z = np.empty((n, n), order="F")
    isuppz = np.empty(2 * n, dtype=_int_dtype)
    m, info = _int(0), _int(0)

    def dsyevr(work, lwork, iwork, liwork):
        _blas.dsyevr_(b"V", b"A", b"L", _ref(_int(n)), _pointer(a), _ref(_int(n)),
                      _ref(_dbl(0.0)), _ref(_dbl(0.0)), _ref(_int(0)), _ref(_int(0)),
                      _ref(_dbl(0.0)), _ref(m), _pointer(values), _pointer(z), _ref(_int(n)),
                      _pointer(isuppz), work, _ref(_int(lwork)), iwork, _ref(_int(liwork)),
                      _ref(info), ctypes.c_size_t(1), ctypes.c_size_t(1), ctypes.c_size_t(1))
        if info.value != 0:
            raise ValueError(f"error code {info.value} from Lapack routine 'dsyevr'")

    size, isize = _dbl(0.0), _int(0)
    dsyevr(_ref(size), -1, _ref(isize), -1)
    work = np.empty(int(size.value))
    iwork = np.empty(isize.value, dtype=_int_dtype)
    dsyevr(_pointer(work), len(work), _pointer(iwork), len(iwork))
    return {"values": values[::-1].copy(), "vectors": np.asfortranarray(z[:, ::-1])}


# ----- Sums, means, variances -----

@njit(cache=True)
def _sum(x):
    total = 0.0
    for v in x:
        total += v
    return total


def r_sum(x):
    """`sum(x)`: the elements added in order. Returns: float."""
    return float(_sum(np.asarray(x, dtype=np.float64).ravel()))


@njit(cache=True)
def _mean(x):
    n = len(x)
    s = 0.0
    for v in x:
        s += v
    s /= n
    if np.isfinite(s):
        t = 0.0
        for v in x:
            t += v - s
        s += t / n
    return s


def r_mean(x, na_rm=False):
    """`mean(x, na.rm = na_rm)` of doubles: the mean with R's second pass.

    Parameters: x (vector), na_rm (drop missing values). Returns: float.
    """
    x = np.asarray(x, dtype=np.float64)
    if na_rm:
        x = x[~np.isnan(x)]
    if len(x) == 0:
        return float("nan")
    return float(_mean(x))


def r_mean_integer(x, na_rm=False):
    """`mean(x, na.rm = na_rm)` of integers or logicals: the sum divided by
    the length, without the second pass R makes for doubles.

    Parameters: x (vector of whole numbers; NaN is a missing value), na_rm.
    Returns: float.
    """
    x = np.asarray(x, dtype=np.float64)
    if na_rm:
        x = x[~np.isnan(x)]
    if len(x) == 0:
        return float("nan")
    return float(_sum(x) / len(x))


@njit(cache=True)
def _cov(X):
    """Covariance matrix of the columns of X, missing values propagated:
    cov_na_1() of R's cov.c.

    R's compiled loop adds the products of deviations eight rows at a time,
    each product rounded before it is added; the rows left over at the end
    (fewer than eight) are added with one rounding each (muladd) on arm64."""
    n, p = X.shape
    blocks = n - n % 8
    out = np.empty((p, p))
    means = np.empty(p)
    missing = np.zeros(p, dtype=np.bool_)
    for i in range(p):
        for k in range(n):
            if np.isnan(X[k, i]):
                missing[i] = True
                break
        if not missing[i]:
            means[i] = _mean(X[:, i])
    for i in range(p):
        for j in range(i + 1):
            if missing[i] or missing[j] or n <= 1:
                out[i, j] = out[j, i] = np.nan
            else:
                total = 0.0
                for k in range(n):
                    if k < blocks:
                        total += (X[k, i] - means[i]) * (X[k, j] - means[j])
                    else:
                        total = muladd(X[k, i] - means[i], X[k, j] - means[j], total)
                out[i, j] = out[j, i] = total / (n - 1)
    return out


def r_var(x, na_rm=False):
    """`var(x, na.rm = na_rm)`: the variance of a vector, or the covariance
    matrix of the columns of a matrix.

    Parameters: x (vector or matrix), na_rm (vector only: drop missing
    values). Returns: float, or Fortran-ordered matrix.
    """
    x = np.asarray(x, dtype=np.float64)
    if x.ndim == 2:
        return np.asfortranarray(_cov(np.asfortranarray(x)))
    if na_rm:
        x = x[~np.isnan(x)]
    return float(_cov(x.reshape(-1, 1))[0, 0])


def r_sd(x, na_rm=False):
    """`sd(x, na.rm = na_rm)`. Returns: float."""
    variance = r_var(x, na_rm)
    return float(np.sqrt(variance)) if variance >= 0 or math.isnan(variance) else float("nan")


@njit(cache=True)
def _col_means(A):
    out = np.empty(A.shape[1])
    for j in range(A.shape[1]):
        total = 0.0
        for i in range(A.shape[0]):
            total += A[i, j]
        out[j] = total / A.shape[0]
    return out


def col_means(A):
    """`colMeans(A)`: each column added in order, divided by its length.

    Parameters: A (matrix). Returns: float array.
    """
    return _col_means(np.asarray(A, dtype=np.float64))


def quantile(x, prob):
    """`quantile(x, probs = prob)[[1]]`, type 7. Returns: float."""
    x = np.sort(np.asarray(x, dtype=np.float64))
    index = 1 + max(len(x) - 1, 0) * prob
    lo, hi = int(np.floor(index)), int(np.ceil(index))
    q = x[lo - 1]
    if index > lo and x[hi - 1] != q:
        h = index - lo
        q = (1 - h) * q + h * x[hi - 1]
    return float(q)


# ----- Smallest and largest, missing values as in R -----

def r_min(*values):
    """`min(...)` of numbers: NaN when one is missing. Returns: float."""
    values = [float(v) for v in values]
    return float("nan") if any(math.isnan(v) for v in values) else min(values)


def r_max(*values):
    """`max(...)` of numbers: NaN when one is missing. Returns: float."""
    values = [float(v) for v in values]
    return float("nan") if any(math.isnan(v) for v in values) else max(values)


def pmin(bound, x):
    """`pmin(bound, x)` for a number and a vector: the smaller of the two
    at every element, NaN where x is missing. Returns: float array."""
    x = np.asarray(x, dtype=np.float64)
    return np.where(np.isnan(x), np.nan, np.where(x > bound, bound, x))


def pmax(bound, x):
    """`pmax(bound, x)`: the larger of the two at every element, NaN where
    one is missing. Parameters: bound (number or array), x (array).
    Returns: float array."""
    x = np.asarray(x, dtype=np.float64)
    bound = np.broadcast_to(np.asarray(bound, dtype=np.float64), x.shape)
    return np.where(np.isnan(x) | np.isnan(bound), np.nan, np.where(x < bound, bound, x))


def tabulate(x, nbins):
    """`tabulate(x, nbins)`: how often each of 1, ..., nbins occurs in x;
    other values are ignored.

    Parameters: x (whole numbers), nbins. Returns: float array.
    """
    x = np.asarray(x, dtype=np.float64)
    x = x[~np.isnan(x)].astype(np.int64)
    x = x[(x >= 1) & (x <= nbins)]
    return np.bincount(x, minlength=int(nbins) + 1)[1:].astype(np.float64)


def count_frequencies_from_zero(x, largest):
    """`c(sum(x == 0), tabulate(x, largest))`: how often each of 0, ...,
    largest occurs in x.

    Parameters: x (whole numbers), largest. Returns: float array.
    """
    x = np.asarray(x, dtype=np.float64)
    return np.concatenate(([float(np.sum(x == 0))], tabulate(x, largest)))


def r_order(x):
    """`order(x)`, zero-based: the stable order of a vector, missing values
    last. Returns: int array."""
    return np.argsort(np.asarray(x, dtype=np.float64), kind="stable")


# ----- Tables in files -----

def read_csv(file, **options):
    """`read.csv(file)`: a csv file as a data frame; an empty cell or NA is
    a missing value.

    Parameters: file (path), options (passed to pandas).
    Returns: pandas DataFrame.
    """
    return pd.read_csv(file, float_precision="round_trip", keep_default_na=False,
                       na_values=["", "NA"], **options)


_EXACT_POWER = 22  # largest power of ten that is an exact double


def _scaled_integer(r, kp, digits):
    """r / 10^kp rounded to the nearest integer (ties to even), and kp, one
    lower when the quotient has fewer than `digits` digits: a step of R's
    scientific()."""
    if abs(kp) <= _EXACT_POWER:
        scaled = r / 10.0 ** kp if kp >= 0 else r * 10.0 ** -kp
    elif kp <= -308:
        scaled = (r * 1e+303) / math.pow(10.0, kp + 303)
    else:
        scaled = r / math.pow(10.0, kp)
    if scaled < 10.0 ** (digits - 1):
        scaled *= 10.0
        kp -= 1
    return float(round(scaled)), kp


def format_number(x, digits=15):
    """A number as R writes it with `digits` significant digits
    (`as.character()`, `write.csv()`).

    R finds the significant digits a number needs by scaling it to an
    integer of `digits` digits and dropping the zeros at its end, then
    writes it in fixed or scientific notation, whichever is shorter.

    Parameters: x (float), digits (significant digits, at most 15).
    Returns: the text; NA for a missing value.
    """
    x = float(x)
    if math.isnan(x):
        return "NA"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    if x == 0:
        return "0"
    negative, r = int(x < 0), abs(x)
    kp = math.floor(math.log10(r)) - digits + 1
    alpha, kp = _scaled_integer(r, kp, digits)
    significant = digits
    for _ in range(digits):
        alpha /= 10.0
        if alpha != math.floor(alpha):
            break
        significant -= 1
    if significant == 0:
        significant = 1
        kp += 1
    power = kp + digits - 1
    left = power + 1
    cut = min(max(digits - power, 0), _EXACT_POWER)  # the decimal place rounding cuts off
    if 0 < power <= _EXACT_POWER and r < 10.0 ** power - 0.5 / 10.0 ** cut:
        left -= 1
    decimals = max(0, significant - power - 1)
    fixed_width = negative + max(left, 1) + decimals + (decimals > 0)
    mantissa = significant - 1
    scientific_width = negative + (mantissa > 0) + mantissa + 4 + (2 if abs(power) >= 100 else 1)
    return f"{x:.{decimals}f}" if fixed_width <= scientific_width else f"{x:.{mantissa}e}"


def write_csv(table, file):
    """`write.csv(table, file, row.names = FALSE)`: numbers with 15
    significant digits as R writes them, text in quotes, NA for a missing
    value.

    Parameters: table (DataFrame), file (path). Returns: None.
    """
    def quoted(text):
        return '"' + str(text).replace('"', '""') + '"'

    columns = []
    for name in table.columns:
        values = table[name]
        if pd.api.types.is_bool_dtype(values):
            columns.append(["TRUE" if value else "FALSE" for value in values])
        elif pd.api.types.is_integer_dtype(values):
            columns.append([str(value) for value in values])
        elif pd.api.types.is_float_dtype(values):
            columns.append([format_number(float(value)) for value in values])
        else:
            columns.append(["NA" if value is None or (isinstance(value, float) and math.isnan(value))
                            else quoted(value) for value in values])
    with open(file, "w") as f:
        f.write(",".join(quoted(name) for name in table.columns) + "\n")
        for row in zip(*columns):
            f.write(",".join(row) + "\n")


# ----- Parallel processes -----

_job = None


def _run_job(item):
    try:
        return _job(item)
    except Exception as error:  # as mclapply: the error is the result
        return error


def parallel_map(function, items, cores=1):
    """`parallel::mclapply(items, function, mc.cores = cores,
    mc.preschedule = FALSE)`.

    The processes are forked, so each starts with the memory and the random
    number stream of this one. A job that fails gives its exception as its
    result, as mclapply gives the error.

    Parameters: function (of one item), items (list), cores (number of
    parallel processes).
    Returns: list of the results, in the order of `items`.
    """
    global _job
    items = list(items)
    cores = int(cores)
    _job = function
    if cores <= 1 or len(items) <= 1:
        return [_run_job(item) for item in items]
    with multiprocessing.get_context("fork").Pool(min(cores, len(items))) as pool:
        return pool.map(_run_job, items, chunksize=1)
