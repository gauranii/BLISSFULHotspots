#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# File:    checks/run_checks.sh
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Build every table and figure with the default settings of its
#          script, compare each with the paper by code, and write the
#          comparison as Markdown.
# Usage:   bash checks/run_checks.sh [RESULTS.md]
#          Run from the repository root. Outputs go to results/ and figures/
#          (or under LOCALFDR_OUTPUT). The comparison is written to the file
#          given, checks/RESULTS.md by default. Exit status 1 if a check fails.
#          CORES sets the parallel processes of Table 3 (default 4).
# -----------------------------------------------------------------------------
set -euo pipefail
python=${PYTHON:-python3}  # a Python 3 interpreter (standard library only)
output=${LOCALFDR_OUTPUT:-$(pwd)}
report=${1:-checks/RESULTS.md}

for script in scripts/*.R; do
  case $script in
    *table3*) Rscript "$script" --cores="${CORES:-4}" > /dev/null ;;
    *) Rscript "$script" > /dev/null ;;
  esac
done

status=0
"$python" checks/check_tables.py "$output/results" "$output/results/checks" || status=1
"$python" checks/check_figures.py archive/figures "$output/figures" "$output/results" \
  "$output/results/checks" || status=1
"$python" checks/write_results.py "$output/results/checks" > "$report"
exit $status
