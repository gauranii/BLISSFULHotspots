# -----------------------------------------------------------------------------
# File:    checks/write_rcpp.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Write the comparison of the R and the Rcpp version as Markdown
#          (checks/RCPP.md) from the tables of checks/check_python.py.
# Usage:   python3 checks/write_rcpp.py <rcpp_archive.csv> <rcpp_rerun.csv>
#          [seconds to compile] > checks/RCPP.md
#          A table that does not exist is left out.
# -----------------------------------------------------------------------------
import csv
import os
import sys

NAMES = {"table1": "Table 1", "table2": "Table 2", "table3": "Table 3", "fig1": "Figure 1",
         "figS1": "Figure S1", "figS2": "Figure S2", "figS3": "Figure S3"}
INTRO = {
    "archive": """## With the defaults of the scripts

Both versions run with the defaults of their scripts, as `nix flake check`
runs them: Tables 1 and 2 and Figures S1 and S3 from the archived runs,
Table 3 from the data (18 seeded chains of 11000 iterations), Figure 1 from
the data and Figure S2 from seeded simulated data.
""",
    "rerun": """## Recomputed by simulation

Both versions run Tables 1 and 2 with `--rerun` for the first replications of
every setting and method, with the seeds of the full runs
(`bash checks/check_rcpp.sh rerun <work dir>`). Compared are the tables and
the stored file of every replication (rejections, false and true rejections
of each data set and method).
""",
}


def section(file):
    rows = list(csv.DictReader(open(file)))
    mode = rows[0]["mode"]
    identical = sum(row["verdict"] == "identical" for row in rows)
    lines = [INTRO[mode], f"{identical} of {len(rows)} outputs identical in every value "
             f"({rows[0]['cores']} cores).", "",
             "| Output | From | Files | Values | Not identical | Largest difference | Verdict | "
             "R (s) | Rcpp (s) |",
             "|---|---|---|---|---|---|---|---|---|"]
    for row in rows:
        lines.append(f"| {NAMES[row['output']]} | {row['source']} | {row['files']} | "
                     f"{row['values']} | {row['not_identical']} | {row['max_abs_diff']} | "
                     f"{row['verdict']} | {row['r_seconds']} | {row['rcpp_seconds']} |")
    problems = [f"- {NAMES[row['output']]}: {row['problems']}" for row in rows if row["problems"]]
    if problems:
        lines += ["", "Differences:", ""] + problems
    return "\n".join(lines) + "\n"


compiling = (f"\nCompiling the C++ files took {sys.argv[3]} s, once; the Rcpp times are those of\n"
             "the scripts with the compiled files in place."
             if len(sys.argv) > 3 and sys.argv[3] else "")
print(f"""# R and Rcpp

Written by `checks/check_rcpp.sh`: the outputs of the R scripts (`scripts/`)
and of the Rcpp scripts (`rcpp/scripts/`) compared by code
(`checks/check_python.py`, the comparison also used for the Python version).
Every csv file under `results/` is compared cell by cell, without tolerance:
the same header, the same rows, every cell the same text or the same number.
Every figure is compared by the marks read from the two PDF files
(`checks/figure_compare.py`): the same number of bars, lines or tiles in the
same order, with the same fills and positions. "Values" counts the compared
cells and mark coordinates. Times are wall-clock seconds on an 8 core Apple
M2 laptop shared with other jobs.{compiling}

The Rcpp scripts run the R scripts with the C++ functions of `rcpp/src/` in
place of the R functions of the same names: the distributions, the whole
Gibbs sampler and the simulated data. The figures themselves, the empirical
Bayes methods and the tables are made by the same R code in both versions.

The tests in `rcpp/tests/` (also run by `nix flake check`) call every C++
function and its R function with the same arguments under the same seed and
expect identical values and the same state of the random number generator
afterwards: the ZIGP probabilities and draws, the Dirichlet and multinomial
draws, every step of the sampler, the adapted proposal, the simulated data,
and seeded Gibbs chains for the three versions of the sampler, the four
alternatives, both switches and three data sets.
""")
for file in sys.argv[1:3]:
    if os.path.exists(file):
        print(section(file))
