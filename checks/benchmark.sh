#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# File:    checks/benchmark.sh
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Time the R code against the Rcpp code on the same seeded work and
#          check that both give the same outputs:
#          chain   one Gibbs chain at the settings of Table 3 (domain cd00063,
#                  the largest, NP alternative), 1 core;
#          table3  the whole of Table 3 (18 chains of 11000 iterations),
#                  1 and 2 cores;
#          table1  one replication of every setting and method of Table 1
#                  (--rerun --runs=1), 1 and 2 cores;
#          each RUNS times (default 5), R and Rcpp in turn. Then, once and on
#          1 core, one replication of every simulation study of the paper
#          (Tables 1 and 2, Figures S1 and S3, in this order, replications
#          shared by two studies run once), from which the time of the full
#          runs of 1000 replications is estimated.
# Outputs: <work dir>/times.csv (one line per timed run) and
#          <work dir>/identical.csv (one line per comparison of an R output
#          with the Rcpp output of the same run); checks/write_benchmark.R
#          turns them into checks/BENCHMARK.md.
# Usage:   bash checks/benchmark.sh <work dir>
#          Run from the repository root, inside `nix develop`, on a machine
#          with no other jobs. RUNS sets the runs per cell. Finished runs are
#          kept, so an interrupted benchmark continues.
#          Exit status 1 if an R output and its Rcpp output differ.
# -----------------------------------------------------------------------------
set -euo pipefail
work=$1
runs=${RUNS:-5}
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export LOCALFDR_RCPP_CACHE=$work/compiled
mkdir -p "$work"
times=$work/times.csv
same=$work/identical.csv
[ -e "$times" ] || echo "cell,impl,cores,run,seconds" > "$times"
[ -e "$same" ] || echo "cell,cores,run,files,identical" > "$same"

# Compile the C++ files once, so that no time below includes it.
if [ ! -e "$work/compiled.done" ]; then
  Rscript -e 'source("R/setup.R"); source("rcpp/R/load.R")' > "$work/compiled.log" 2>&1
  touch "$work/compiled.done"
fi

script() {  # script <impl> <name>
  if [ "$1" = r ]; then echo "scripts/$2.R"; else echo "rcpp/scripts/$2.R"; fi
}

# Run one cell once with one implementation, unless already done.
run_cell() {  # run_cell <cell> <impl> <cores> <run>
  local cell=$1 impl=$2 cores=$3 run=$4
  local dir=$work/$cell/$impl-c$cores-run$run
  [ -e "$dir.done" ] && return 0
  rm -rf "$dir"; mkdir -p "$dir"
  local start=$EPOCHREALTIME seconds
  case $cell in
    chain)
      seconds=$(Rscript checks/benchmark_chain.R "$impl" cd00063 NP "$dir/fdr.csv" | cut -d, -f4) ;;
    table3)
      LOCALFDR_OUTPUT=$dir Rscript "$(script "$impl" table3_protein_domains)" --cores="$cores" \
        > "$dir.log" 2>&1 ;;
    table1)
      LOCALFDR_OUTPUT=$dir Rscript "$(script "$impl" table1_simulation)" \
        --rerun --runs=1 --replications=1 --cores="$cores" > "$dir.log" 2>&1 ;;
  esac
  [ -n "${seconds:-}" ] || seconds=$(awk -v a="$start" -v b="$EPOCHREALTIME" 'BEGIN { printf "%.3f", b - a }')
  echo "$cell,$impl,$cores,$run,$seconds" >> "$times"
  touch "$dir.done"
  echo "$cell $impl cores=$cores run=$run: $seconds s"
}

# Compare the R output of a run with its Rcpp output, file by file and byte by byte.
compare() {  # compare <cell> <cores> <run>
  local a=$work/$1/r-c$2-run$3 b=$work/$1/rcpp-c$2-run$3 files ok=yes
  grep -q "^$1,$2,$3," "$same" && return 0
  files=$( (cd "$a" && find . -type f | sort) )
  [ "$files" = "$( (cd "$b" && find . -type f | sort) )" ] || ok=no
  diff -r -q "$a" "$b" > /dev/null || ok=no
  echo "$1,$2,$3,$(echo "$files" | wc -l | tr -d ' '),$ok" >> "$same"
}

for run in $(seq 1 "$runs"); do
  for impl in r rcpp; do run_cell chain "$impl" 1 "$run"; done
  compare chain 1 "$run"
  for cores in 1 2; do
    for cell in table3 table1; do
      for impl in r rcpp; do run_cell "$cell" "$impl" "$cores" "$run"; done
      compare "$cell" "$cores" "$run"
    done
  done
done

# One replication of every simulation study, 1 core, one folder per
# implementation so that shared replications are run once, as in a full run.
for impl in r rcpp; do
  for name in table1_simulation table2_prior_sensitivity figS1_fdr_tpr_well_separated figS3_fdr_table; do
    dir=$work/study/$impl
    [ -e "$dir/$name.done" ] && continue
    mkdir -p "$dir"
    start=$EPOCHREALTIME
    LOCALFDR_OUTPUT=$dir Rscript "$(script "$impl" "$name")" \
      --rerun --runs=1 --replications=1 --cores=1 > "$dir/$name.log" 2>&1
    seconds=$(awk -v a="$start" -v b="$EPOCHREALTIME" 'BEGIN { printf "%.3f", b - a }')
    echo "study-$name,$impl,1,1,$seconds" >> "$times"
    touch "$dir/$name.done"
    echo "study $name $impl: $seconds s"
  done
done
if ! grep -q "^study," "$same"; then
  files=$(cd "$work/study/r/results" && find . -type f | wc -l | tr -d ' ')
  ok=yes; diff -r -q "$work/study/r/results" "$work/study/rcpp/results" > /dev/null || ok=no
  echo "study,1,1,$files,$ok" >> "$same"
fi

! grep -q ",no$" "$same"
