# -----------------------------------------------------------------------------
# File:    checks/write_python.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Write the comparison of the R and the Python implementation as
#          Markdown (checks/PYTHON.md) from the tables of
#          checks/check_python.py.
# Usage:   python3 checks/write_python.py <python_archive.csv>
#          <python_rerun.csv> > checks/PYTHON.md
#          A table that does not exist is left out.
# -----------------------------------------------------------------------------
import csv
import os
import sys

NAMES = {"table1": "Table 1", "table2": "Table 2", "table3": "Table 3", "fig1": "Figure 1",
         "figS1": "Figure S1", "figS2": "Figure S2", "figS3": "Figure S3"}
INTRO = {
    "archive": """## With the defaults of the scripts

Both implementations run with the defaults of their scripts, as
`nix flake check` runs them: Tables 1 and 2 and Figures S1 and S3 from the
archived runs, Table 3 from the data (18 seeded chains of 11000 iterations),
Figure 1 from the data and Figure S2 from seeded simulated data.
""",
    "rerun": """## Recomputed by simulation

Both implementations run Tables 1 and 2 with `--rerun` for the first
replications of every setting and method, with the seeds of the full runs
(`bash checks/check_python.sh rerun <work dir>`). Compared are the tables and
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
             "R (s) | Python (s) |",
             "|---|---|---|---|---|---|---|---|---|"]
    for row in rows:
        lines.append(f"| {NAMES[row['output']]} | {row['source']} | {row['files']} | "
                     f"{row['values']} | {row['not_identical']} | {row['max_abs_diff']} | "
                     f"{row['verdict']} | {row['r_seconds']} | {row['python_seconds']} |")
    problems = [f"- {NAMES[row['output']]}: {row['problems']}" for row in rows if row["problems"]]
    if problems:
        lines += ["", "Differences:", ""] + problems
    return "\n".join(lines) + "\n"


print("""# R and Python

Written by `checks/check_python.sh`: the outputs of the R scripts (`scripts/`)
and of the Python scripts (`python/scripts/`) compared by code
(`checks/check_python.py`). Every csv file under `results/` is compared cell
by cell, without tolerance: the same header, the same rows, every cell the
same text or the same number. Every figure is compared by the marks read from
the two PDF files (`checks/figure_compare.py`): the same number of bars, lines
or tiles in the same order and with the same fills; bar heights, bar positions
and line vertices agree to 0.02 pt of the R page after one straight-line map
per axis, since R and matplotlib lay their pages out differently and R prints
its coordinates with two decimals. "Values" counts the compared cells and
mark coordinates. Times are wall-clock seconds on an 8 core Apple M2 laptop
shared with other jobs; the Python times include loading its compiled (numba)
functions.

The tests in `python/tests/` (also run by `nix flake check`) compare with R,
bit for bit: the random numbers of `python/localfdr/rng.py` (uniform, normal,
exponential, gamma, beta, binomial, Poisson, geometric, `sample.int`) for
several seeds and over every branch of the algorithms; `lgamma`, the logistic
functions and `round`; the Cholesky factor, the eigenvalues and the
covariances behind the adapted proposal; the ZIGP probabilities and draws;
the empirical Bayes fit of five simulated data sets and of the six domains;
and every stored draw of seeded Gibbs chains for the three versions of the
sampler, the four alternatives and three data sets.
""")
for file in sys.argv[1:3]:
    if os.path.exists(file):
        print(section(file))
