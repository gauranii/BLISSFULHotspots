# -----------------------------------------------------------------------------
# File:    python/tests/test_rbase.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Tests of rbase.py: every function on small inputs with a known
#          answer, and the same doubles and the same text as R for the
#          decompositions, the variances and the numbers written to csv.
# Usage:   python -m pytest python/tests/test_rbase.py
# -----------------------------------------------------------------------------
import math

import numpy as np
import pandas as pd
import pytest

from localfdr import rbase


def test_as_matrix_and_matprod():
    assert rbase.as_matrix([1, 2, 3]).shape == (3, 1)
    A = np.array([[1.0, 2.0], [3.0, 4.0]])
    np.testing.assert_array_equal(rbase.matprod(A, A), [[7.0, 10.0], [15.0, 22.0]])
    np.testing.assert_array_equal(rbase.matprod(A, [1.0, 1.0]), [[3.0], [7.0]])
    np.testing.assert_array_equal(rbase.matprod(np.array([[1.0, 1.0]]), A), [[4.0, 6.0]])
    with_missing = rbase.matprod(np.array([[1.0, np.nan], [0.0, 1.0]]), A)
    assert np.isnan(with_missing[0]).all() and list(with_missing[1]) == [3.0, 4.0]


def test_chol():
    A = np.array([[4.0, 2.0], [2.0, 3.0]])
    R = rbase.chol(A)
    np.testing.assert_allclose(R, [[2.0, 1.0], [0.0, math.sqrt(2.0)]], rtol=1e-15)
    np.testing.assert_allclose(R.T @ R, A, rtol=1e-15)
    with pytest.raises(ValueError):
        rbase.chol(np.array([[1.0, 2.0], [2.0, 1.0]]))
    with pytest.raises(ValueError):
        rbase.chol(np.ones((2, 3)))


def test_eigen_symmetric():
    e = rbase.eigen_symmetric(np.diag([1.0, 3.0, 2.0]))
    np.testing.assert_allclose(e["values"], [3.0, 2.0, 1.0], rtol=1e-15)
    np.testing.assert_allclose(np.abs(e["vectors"]), [[0, 0, 1], [1, 0, 0], [0, 1, 0]], atol=1e-15)
    S = np.array([[2.0, 1.0], [1.0, 2.0]])
    e = rbase.eigen_symmetric(S)
    np.testing.assert_allclose(e["values"], [3.0, 1.0], rtol=1e-14)
    np.testing.assert_allclose(e["vectors"] @ np.diag(e["values"]) @ e["vectors"].T, S, atol=1e-14)
    with pytest.raises(ValueError):
        rbase.eigen_symmetric(np.array([[1.0, np.nan], [np.nan, 1.0]]))
    with pytest.raises(ValueError):
        rbase.eigen_symmetric(np.array([[1.0, 2.0], [0.0, 1.0]]))


def test_sums_and_means():
    assert rbase.r_sum([1, 2, 3.5]) == 6.5
    assert rbase.r_sum([]) == 0.0
    assert rbase.r_mean([1.0, 2.0, 6.0]) == 3.0
    assert rbase.r_mean([1.0, np.nan, 3.0], na_rm=True) == 2.0
    assert math.isnan(rbase.r_mean([1.0, np.nan])) and math.isnan(rbase.r_mean([]))
    assert rbase.r_mean_integer([1, 2, 4]) == 7 / 3
    assert rbase.r_mean_integer([True, False, False, True]) == 0.5
    assert rbase.r_mean_integer([1, np.nan, 3], na_rm=True) == 2.0
    assert math.isnan(rbase.r_mean_integer([], na_rm=True))
    np.testing.assert_array_equal(rbase.col_means(np.array([[1.0, 2.0], [3.0, 6.0]])), [2.0, 4.0])


def test_variances():
    assert rbase.r_var([1.0, 2.0, 3.0, 4.0]) == pytest.approx(5 / 3, rel=1e-15)
    assert rbase.r_sd([2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]) == pytest.approx(math.sqrt(32 / 7))
    assert math.isnan(rbase.r_var([1.0])) and math.isnan(rbase.r_sd([1.0]))
    assert math.isnan(rbase.r_var([1.0, np.nan, 2.0]))
    assert rbase.r_var([1.0, np.nan, 3.0], na_rm=True) == 2.0
    X = np.array([[1.0, 2.0], [2.0, 4.0], [3.0, 7.0], [4.0, 7.0]])
    np.testing.assert_allclose(rbase.r_var(X), np.cov(X, rowvar=False), rtol=1e-15)
    X[1, 1] = np.nan
    V = rbase.r_var(X)
    assert V[0, 0] == pytest.approx(5 / 3) and np.isnan(V[0, 1]) and np.isnan(V[1, 1])


def test_quantile():
    assert rbase.quantile([1.0, 2.0, 3.0, 4.0], 0.75) == 3.25
    assert rbase.quantile([5.0, 1.0, 3.0], 0.5) == 3.0
    assert rbase.quantile([2.0], 0.75) == 2.0


def test_smallest_and_largest_with_missing_values():
    assert rbase.r_min(3, 1.5) == 1.5 and rbase.r_max(3, 1.5) == 3.0
    assert math.isnan(rbase.r_min(1, math.nan)) and math.isnan(rbase.r_max(math.nan, 1))
    out = rbase.pmin(1.0, [0.5, 2.0, np.nan, np.inf])
    np.testing.assert_array_equal(out[[0, 1, 3]], [0.5, 1.0, 1.0])
    assert np.isnan(out[2])
    out = rbase.pmax([1.0, 0.0, 1.0], [0.5, 0.2, np.nan])
    np.testing.assert_array_equal(out[:2], [1.0, 0.2])
    assert np.isnan(out[2])


def test_tabulate_and_order():
    np.testing.assert_array_equal(rbase.tabulate([1, 3, 3, 0, 7], 4), [1, 0, 2, 0])
    np.testing.assert_array_equal(rbase.count_frequencies_from_zero([0, 0, 2, 3, 3, 3], 3),
                                  [2, 0, 1, 3])
    np.testing.assert_array_equal(rbase.r_order([0.3, 0.1, 0.3, 0.0]), [3, 1, 0, 2])
    np.testing.assert_array_equal(rbase.r_order([1.0, np.nan, 0.5]), [2, 0, 1])


@pytest.mark.parametrize("value, text", [
    (0.0001, "1e-04"), (100000.0, "1e+05"), (123456.0, "123456"), (0.1 + 0.2, "0.3"),
    (1 / 3, "0.333333333333333"), (0.0432, "0.0432"), (3e-04, "3e-04"), (0.0011, "0.0011"),
    (-2.5, "-2.5"), (50.0, "50"), (1e-10, "1e-10"), (1234567.891, "1234567.891"),
    (float("nan"), "NA"), (float("inf"), "Inf"), (0.0, "0"), (1e15, "1e+15"), (0.02, "0.02")])
def test_format_number(value, text):
    assert rbase.format_number(value) == text


def test_csv_files(tmp_path):
    table = pd.DataFrame({"name": ["a", 'b "c"', None], "count": [1, 2, 3],
                          "value": [0.5, np.nan, 1e-04], "flag": [True, False, True]})
    file = tmp_path / "table.csv"
    rbase.write_csv(table, file)
    assert file.read_text().splitlines() == [
        '"name","count","value","flag"', '"a",1,0.5,TRUE', '"b ""c""",2,NA,FALSE', 'NA,3,1e-04,TRUE']
    back = rbase.read_csv(file)
    assert list(back["count"]) == [1, 2, 3] and math.isnan(back["value"][1])
    assert back["value"][2] == 1e-04 and pd.isna(back["name"][2])


def _square(item):
    if item == 3:
        raise ValueError("three")
    return item * item


def test_parallel_map():
    assert rbase.parallel_map(_square, [1, 2, 4], cores=1) == [1, 4, 16]
    assert rbase.parallel_map(_square, [1, 2, 4, 5], cores=2) == [1, 4, 16, 25]
    results = rbase.parallel_map(_square, [2, 3], cores=2)
    assert results[0] == 4 and isinstance(results[1], ValueError)


@pytest.mark.against_r
def test_same_doubles_and_text_as_r(run_r, reader, tmp_path):
    numbers = tmp_path / "numbers.csv"
    r = reader(run_r(r"""
set.seed(9); H <- matrix(rnorm(300), 100, 3) %%*%% matrix(c(1, .5, .2, 0, 1, .3, 0, 0, .1), 3)
H[60:75, ] <- H[rep(59, 16), ]
put(H); V <- var(H[50:100, ]); put(V); s <- 50 * V / 51
e <- eigen(s); put(e$values); put(e$vectors); P <- positive_definite(s); put(P); put(chol(P))
put(e$vectors %%*%% diag(c(1e-3, 2e-3, 3e-3), 3)); put(t(chol(P)) %%*%% c(0.3, -1.1, 0.7))
A <- H[1:3, ]; A[2, 1] <- Inf; put(A %%*%% H[4:6, ])
for (n in 2:40) { x <- rnorm(n) * 10^(n %%%% 5); put(x); put(c(var(x), sd(x), mean(x), quantile(x, 0.75), sum(x))) }
put(mean(c(3L, 4L, 4L))); put(colMeans(H))
z <- c(H[1:40, 1], round(H[1:40, 2], 4), 1 / 3, 1e5, 1e-4, 123456, 1e15, 1e-15, 0.1 + 0.2, 100, 3e-4,
       2^60, 1e22, 1e23, 5e-324, -1.5e-7, 0.99999999999999989, 1e5 + 0.1, 123456789012345678)
put(z); write.csv(data.frame(z = z, text = "a"), "%s", row.names = FALSE)
""" % numbers))
    H = r.take(300).reshape(100, 3, order="F")
    V = rbase.r_var(H[49:100])
    r.same("var of a matrix", V)
    s = 50 * V / 51
    e = rbase.eigen_symmetric(s)
    r.same("eigenvalues", e["values"])
    r.same("eigenvectors", e["vectors"])
    from localfdr.bayes_fdr import positive_definite
    P = positive_definite(s)
    r.same("positive definite", P)
    r.same("chol", rbase.chol(P))
    r.same("dgemm", rbase.matprod(e["vectors"], np.diag([1e-3, 2e-3, 3e-3])))
    r.same("dgemv", rbase.matprod(rbase.chol(P).T, [0.3, -1.1, 0.7]))
    A = H[0:3].copy()
    A[1, 0] = np.inf
    r.same("product with an infinite value", rbase.matprod(A, H[3:6]))
    for n in range(2, 41):
        x = r.take(n)
        r.same(f"var, sd, mean, quantile, sum of {n}",
               [rbase.r_var(x), rbase.r_sd(x), rbase.r_mean(x), rbase.quantile(x, 0.75),
                rbase.r_sum(x)])
    r.same("mean of integers", [rbase.r_mean_integer([3, 4, 4])])
    r.same("colMeans", rbase.col_means(H))
    z = r.take(97)
    r.done()
    written = numbers.read_text().splitlines()[1:]
    assert [rbase.format_number(value) + ',"a"' for value in z] == written
