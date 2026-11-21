# -----------------------------------------------------------------------------
# File:    checks/check_python.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Compare the outputs of the R scripts with those of the Python
#          scripts, output by output: every cell of every csv file under
#          results/ (tables, plotted values, stored replications) and the
#          marks of every figure (bars, lines, tiles; figure_compare.py).
# Usage:   python3 checks/check_python.py archive|rerun <R output dir>
#          <Python output dir> <comparison csv> [replications]
# Inputs:  <dir>/results/, <dir>/figures/, <dir>/logs/<script>.done (time and
#          cores) of both implementations, written by checks/check_python.sh.
# Outputs: the comparison csv, one row per output of the paper: files, cells
#          and marks compared, how many are not identical, the largest
#          difference, the verdict and the run times. Prints the same; exit
#          status 1 unless every output is identical.
#          Tolerance: none for the csv files (the same header, the same rows,
#          every cell the same text or the same number). The marks of a
#          figure must correspond one to one, with the same fills; their
#          positions agree to 0.02 pt of the R page after one straight-line
#          map per axis (R prints its coordinates with two decimals).
# -----------------------------------------------------------------------------
import csv
import glob
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_compare import (COLOUR_TOLERANCE, POSITION_TOLERANCE, compare_bars,  # noqa: E402
                            compare_lines, compare_tiles)

mode, r_dir, py_dir, out_file = sys.argv[1:5]
replications = sys.argv[5] if len(sys.argv) > 5 else "2"

# Script, csv patterns (under results/), figure patterns (under figures/),
# how the marks of the figures are compared, and where the output comes from.
ARCHIVE = {
    "table1": ("table1_simulation", ["Table1_*.csv"], [], None, "archived runs"),
    "table2": ("table2_prior_sensitivity", ["Table2_*.csv"], [], None, "archived runs"),
    "table3": ("table3_protein_domains", ["Table3_*.csv"], [], None, "computed (seeded chains)"),
    "fig1": ("fig1_domain_barplots", ["Fig1_*.csv"], ["Fig1*.pdf"], "bars", "computed"),
    "figS1": ("figS1_fdr_tpr_well_separated", ["FigS1_*.csv"], ["FigS1*.pdf"], "lines",
              "archived runs"),
    "figS2": ("figS2_histograms", ["FigS2_*.csv"], ["FigS2*.pdf"], "bars",
              "computed (seeded data)"),
    "figS3": ("figS3_fdr_table", ["FigS3_*.csv"], ["FigS3*.pdf"], "tiles", "archived runs"),
}
RERUN = {
    "table1": ("table1_simulation", ["Table1_*.csv", "replications/*C10_*/*gamma1_*/*.csv",
                                     "replications/*C10_*/GP_*/*.csv", "replications/*C10_*/EB_*/*.csv",
                                     "replications/*C5_*/*/*.csv"], [], None,
               f"{replications} replications"),
    "table2": ("table2_prior_sensitivity", ["Table2_*.csv", "replications/*C10_*/NP_*/*.csv",
                                            "replications/*C10_*/SP_*/*.csv"], [], None,
               f"{replications} replications"),
}
OUTPUTS = ARCHIVE if mode == "archive" else RERUN
VERTICES = 11  # points of a line of Figure S1 (values of N)


def files_of(folder, patterns):
    """Relative paths of the files under a folder that match the patterns."""
    found = set()
    for pattern in patterns:
        found.update(glob.glob(os.path.join(folder, pattern)))
    return sorted(os.path.relpath(path, folder) for path in found)


def number(cell):
    """The number of a cell, NaN for a missing value, None for text."""
    if cell in ("NA", "NaN"):
        return math.nan
    try:
        return float(cell)
    except ValueError:
        return None


class Comparison:
    """Running totals of the comparison of one output."""

    def __init__(self):
        self.files = self.values = self.unequal = 0
        self.max_abs = 0.0
        self.problems = []

    def problem(self, text):
        self.problems.append(text)

    def csv_file(self, r_file, py_file, name):
        """Same header, same number of rows, every cell the same."""
        with open(r_file, newline="") as f:
            r_rows = list(csv.reader(f))
        with open(py_file, newline="") as f:
            py_rows = list(csv.reader(f))
        self.files += 1
        if r_rows[0] != py_rows[0] or len(r_rows) != len(py_rows):
            return self.problem(f"{name}: different columns or number of rows")
        for r_row, py_row in zip(r_rows[1:], py_rows[1:]):
            for a, b in zip(r_row, py_row):
                self.values += 1
                x, y = number(a), number(b)
                if x is None or y is None:
                    same = a == b
                else:
                    same = x == y or (math.isnan(x) and math.isnan(y))
                    if not same and not (math.isnan(x) or math.isnan(y)):
                        self.max_abs = max(self.max_abs, abs(x - y))
                if not same:
                    self.unequal += 1
                    self.problem(f"{name}: cells differ")

    def figure(self, r_file, py_file, name, kind):
        """The marks of a figure correspond and agree within the tolerance."""
        compare = {"bars": compare_bars, "tiles": compare_tiles,
                   "lines": lambda a, b: compare_lines(a, b, VERTICES)}[kind]
        in_r, in_py, compared, difference = compare(r_file, py_file)
        self.files += 1
        self.values += compared
        tolerance = COLOUR_TOLERANCE if kind == "tiles" else POSITION_TOLERANCE
        if difference is None:
            self.unequal += max(in_r, in_py, 1)
            self.problem(f"{name}: {in_r} and {in_py} {kind} that do not correspond")
        elif difference > tolerance:
            self.unequal += compared
            self.problem(f"{name}: {kind} differ by {difference:.3g}")


def run_time(folder, script):
    """Seconds and cores of a script, from its .done file."""
    try:
        with open(os.path.join(folder, "logs", script + ".done")) as f:
            return f.read().strip().split(",")
    except OSError:
        return ["", ""]


rows = []
for output, (script, tables, figures, kind, source) in OUTPUTS.items():
    total = Comparison()
    for folder, patterns in (("results", tables), ("figures", figures)):
        r_files = files_of(os.path.join(r_dir, folder), patterns)
        py_files = files_of(os.path.join(py_dir, folder), patterns)
        if r_files != py_files:
            total.problem(f"the two implementations wrote different files in {folder}/")
        for file in sorted(set(r_files) & set(py_files)):
            paths = os.path.join(r_dir, folder, file), os.path.join(py_dir, folder, file)
            if folder == "results":
                total.csv_file(*paths, file)
            else:
                total.figure(*paths, file, kind)
    if total.files == 0:
        total.problem("no files to compare")
    verdict = "DIFFERENT" if total.problems else "identical"
    r_time, py_time = run_time(r_dir, script), run_time(py_dir, script)
    rows.append({"output": output, "mode": mode, "source": source, "files": total.files,
                 "values": total.values, "not_identical": total.unequal,
                 "max_abs_diff": f"{total.max_abs:.3g}", "verdict": verdict,
                 "r_seconds": r_time[0], "python_seconds": py_time[0], "cores": r_time[1],
                 "problems": "; ".join(sorted(set(total.problems))[:5])})
    print(f"{output:7} {mode:8} files {total.files:4d} values {total.values:6d} "
          f"not identical {total.unequal:5d}  {verdict}")
    for problem in sorted(set(total.problems))[:5]:
        print("        ", problem)

os.makedirs(os.path.dirname(os.path.abspath(out_file)), exist_ok=True)
with open(out_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
sys.exit(0 if all(row["verdict"] == "identical" for row in rows) else 1)
