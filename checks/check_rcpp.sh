#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# File:    checks/check_rcpp.sh
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Build every table and figure with the R scripts (scripts/) and
#          with the Rcpp scripts (rcpp/scripts/), time both, and compare the
#          two by code (checks/check_python.py): every cell of every csv
#          file and the marks of every figure.
# Usage:   bash checks/check_rcpp.sh archive|rerun <work dir> [RCPP.md]
#          [folder of the comparison tables]
#          Run from the repository root, inside `nix develop`.
#          archive: every script with its defaults: Tables 1 and 2 and
#                   Figures S1 and S3 from the archived runs, Table 3 and
#                   Figures 1 and S2 computed; this is what `nix flake check`
#                   runs.
#          rerun:   Tables 1 and 2 by simulation, the first REPLICATIONS
#                   (default 2) replications of every setting and method, with
#                   the seeds of the full runs.
#          The outputs go to <work dir>/<mode>/{r,rcpp}; finished outputs are
#          kept, so an interrupted run continues. The C++ files are compiled
#          first, into <work dir>/compiled (or LOCALFDR_RCPP_CACHE), and the
#          time this takes is recorded apart from the run times. The
#          comparison is written to rcpp_<mode>.csv in the folder given last
#          (checks/ by default), and both modes together to the Markdown file
#          given (checks/RCPP.md by default).
#          CORES (default 4) is passed to the scripts as --cores; IMPLS
#          (default "r rcpp") chooses the implementations to run.
#          Exit status 1 if a comparison fails.
# -----------------------------------------------------------------------------
set -euo pipefail
mode=$1
work=$2
report=${3:-checks/RCPP.md}
tables=${4:-checks}  # folder of rcpp_<mode>.csv
cores=${CORES:-4}
impls=${IMPLS:-r rcpp}
replications=${REPLICATIONS:-2}
python=${PYTHON:-python3}  # a Python 3 interpreter (standard library only)
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export LOCALFDR_RCPP_CACHE=${LOCALFDR_RCPP_CACHE:-$work/compiled}

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

# Compile the C++ files once, so that no script's time includes it.
case " $impls " in
  *" rcpp "*)
    if [ ! -e "$work/compiled.done" ]; then
      mkdir -p "$work"
      start=$EPOCHREALTIME
      Rscript -e 'source("R/setup.R"); source("rcpp/R/load.R")' > "$work/compiled.log" 2>&1 \
        || { echo "compiling rcpp/src failed, see $work/compiled.log" >&2; exit 1; }
      awk -v a="$start" -v b="$EPOCHREALTIME" 'BEGIN { printf "%.0f\n", b - a }' > "$work/compiled.done"
      echo "compiled rcpp/src: $(cat "$work/compiled.done") s"
    fi ;;
esac

for name in $outputs; do
  for impl in $impls; do
    if [ "$impl" = r ]; then
      timed r "$name" Rscript "scripts/$name.R" $(options "$name")
    else
      timed rcpp "$name" Rscript "rcpp/scripts/$name.R" $(options "$name")
    fi
  done
done

[ "$impls" = "r rcpp" ] || exit 0
status=0
"$python" checks/check_python.py "$mode" "$work/$mode/r" "$work/$mode/rcpp" \
  "$tables/rcpp_$mode.csv" "$replications" rcpp || status=1
other=$([ "$mode" = archive ] && echo rerun || echo archive)
if [ ! -e "$tables/rcpp_$other.csv" ] && [ -e "checks/rcpp_$other.csv" ]; then
  cp "checks/rcpp_$other.csv" "$tables/"
fi
"$python" checks/write_rcpp.py "$tables/rcpp_archive.csv" "$tables/rcpp_rerun.csv" \
  "$(cat "$work/compiled.done")" > "$report"
exit $status
