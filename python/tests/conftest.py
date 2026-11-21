# -----------------------------------------------------------------------------
# File:    python/tests/conftest.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Shared set-up of the unit tests: the import path, the marker against_r
#          (tests that run R and compare with it bit for bit; skipped when
#          Rscript is not on the PATH) and the helper that runs R code.
# Usage:   python -m pytest python/tests/conftest.py
# -----------------------------------------------------------------------------
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "python"))

R_HEADER = r"""
source("R/setup.R")
out <- file(commandArgs(TRUE)[1], "wb")
put <- function(x) writeBin(as.double(x), out)
"""


def pytest_configure(config):
    """Name the marker."""
    config.addinivalue_line("markers", "against_r: runs R and compares with it bit for bit")


def pytest_collection_modifyitems(config, items):
    """Skip the comparisons with R when R is not there."""
    if shutil.which("Rscript") is None:
        for item in items:
            if "against_r" in item.keywords:
                item.add_marker(pytest.mark.skip(reason="Rscript is not on the PATH"))


@pytest.fixture
def run_r(tmp_path):
    """Function that runs R code from the repository root, with R/setup.R
    loaded, and returns the doubles the code wrote with put()."""
    def run(code):
        script, values = tmp_path / "code.R", tmp_path / "values.bin"
        script.write_text(R_HEADER + code + "\nclose(out)\n")
        subprocess.run(["Rscript", str(script), str(values)], cwd=ROOT, check=True,
                       capture_output=True)
        return np.fromfile(values)
    return run


class Reader:
    """Reads the doubles of R in the order they were written and compares."""

    def __init__(self, values):
        self.values, self.position = values, 0

    def take(self, count):
        part = self.values[self.position:self.position + count]
        self.position += count
        return part

    def same(self, label, values):
        """Assert that the next doubles of R are these, bit for bit (a
        missing value equals a missing value)."""
        values = np.asarray(values, dtype=np.float64).ravel(order="F")
        reference = self.take(len(values))
        assert len(reference) == len(values), label
        equal = (values == reference) | (np.isnan(values) & np.isnan(reference))
        assert equal.all(), f"{label}: first difference at element {int(np.argmin(equal))}"

    def done(self):
        assert self.position == len(self.values)


@pytest.fixture
def reader():
    return Reader


@pytest.fixture(autouse=True)
def close_figures():
    """Close the figures a test opened."""
    yield
    if "matplotlib.pyplot" in sys.modules:
        sys.modules["matplotlib.pyplot"].close("all")
