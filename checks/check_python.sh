#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# File:    checks/check_python.sh
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Build every table and figure with the R scripts (scripts/) and
#          with the Python scripts (python/scripts/), time both, and compare
#          the two by code (checks/check_python.py): every cell of every csv
#          file and the marks of every figure.
# Usage:   bash checks/check_python.sh archive|rerun <work dir> [PYTHON.md]
#          [folder of the comparison tables]
#          Run from the repository root, inside `nix develop`, with
#          LOCALFDR_BLAS naming the BLAS library of R (flake.nix sets it).
#          archive: every script with its defaults: Tables 1 and 2 and
#                   Figures S1 and S3 from the archived runs, Table 3 and
#                   Figures 1 and S2 computed; this is what `nix flake check`
#                   runs.
#          rerun:   Tables 1 and 2 by simulation, the first REPLICATIONS
#                   (default 2) replications of every setting and method, with
#                   the seeds of the full runs.
#          The outputs go to <work dir>/<mode>/{r,python}; finished outputs
#          are kept, so an interrupted run continues. The comparison is
#          written to python_<mode>.csv in the folder given last (checks/ by
#          default), and both modes together to the Markdown file given
#          (checks/PYTHON.md by default).
#          CORES (default 4) is passed to the scripts as --cores; IMPLS
#          (default "r python") chooses the implementations to run.
#          Exit status 1 if a comparison fails.
# -----------------------------------------------------------------------------
set -euo pipefail
mode=$1
work=$2
report=${3:-checks/PYTHON.md}
tables=${4:-checks}  # folder of python_<mode>.csv
cores=${CORES:-4}
impls=${IMPLS:-r python}
replications=${REPLICATIONS:-2}
python=${PYTHON:-python}
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

if [ "$mode" = rerun ]; then
  outputs="table1_simulation table2_prior_sensitivity"
else
  outputs="table1_simulation table2_prior_sensitivity table3_protein_domains
fig1_domain_barplots figS1_fdr_tpr_well_separated figS2_histograms figS3_fdr_table"
fi

# Options of one script in the given mode.
options() {
  if [ "$mode" = rerun ]; then
    echo "--rerun --runs=$replications --replications=$replications --cores=$cores"
  elif [ "$1" = table3_protein_domains ]; then
    echo "--cores=$cores"
  fi
}

# Run a command and record its time.
timed() {
  local impl=$1 name=$2; shift 2
  local dir=$work/$mode/$impl
  mkdir -p "$dir/logs"
  [ -e "$dir/logs/$name.done" ] && return 0
  local start=$EPOCHREALTIME
  LOCALFDR_OUTPUT=$dir "$@" > "$dir/logs/$name.log" 2>&1 \
    || { echo "$impl $name failed, see $dir/logs/$name.log" >&2; exit 1; }
  awk -v a="$start" -v b="$EPOCHREALTIME" -v c="$cores" \
    'BEGIN { printf "%.1f,%d\n", b - a, c }' > "$dir/logs/$name.done"
  echo "$impl $name: $(cat "$dir/logs/$name.done") (s, cores)"
}

for name in $outputs; do
  for impl in $impls; do
    if [ "$impl" = r ]; then
      timed r "$name" Rscript "scripts/$name.R" $(options "$name")
    else
      timed python "$name" "$python" "python/scripts/$name.py" $(options "$name")
    fi
  done
done

[ "$impls" = "r python" ] || exit 0
status=0
"$python" checks/check_python.py "$mode" "$work/$mode/r" "$work/$mode/python" \
  "$tables/python_$mode.csv" "$replications" || status=1
other=$([ "$mode" = archive ] && echo rerun || echo archive)
if [ ! -e "$tables/python_$other.csv" ] && [ -e "checks/python_$other.csv" ]; then
  cp "checks/python_$other.csv" "$tables/"
fi
"$python" checks/write_python.py "$tables/python_archive.csv" "$tables/python_rerun.csv" > "$report"
exit $status
