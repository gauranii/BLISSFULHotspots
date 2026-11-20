#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# File:    checks/versions.sh
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Print the versions of the software the flake provides: the nixpkgs
#          revision of flake.lock, R and its packages, and Python (used by the
#          checks, standard library only).
# Usage:   bash checks/versions.sh            print the versions
#          bash checks/versions.sh --check    compare them with the section
#                                             "Software versions" of README.md
#          Run from the repository root, inside `nix develop`.
# Outputs: prints the versions; with --check exit status 1 if README.md differs.
# -----------------------------------------------------------------------------
set -euo pipefail

versions() {
  python3 - <<'PY'
import json
import platform

locked = json.load(open("flake.lock"))["nodes"]["nixpkgs"]["locked"]
print(f"nixpkgs  {locked['owner']}/{locked['repo']} {locked['rev']}")
print()
print(f"Python   {platform.python_version()}")
PY
  echo
  Rscript -e '
    cat(sprintf("R        %s.%s\n", R.version$major, R.version$minor))
    installed <- installed.packages()
    installed <- installed[is.na(installed[, "Priority"]) | installed[, "Priority"] != "base", ]
    installed <- installed[order(tolower(rownames(installed))), ]
    cat(sprintf("  %-12s %s\n", rownames(installed), installed[, "Version"]), sep = "")'
}

if [ "${1:-}" = --check ]; then
  diff <(versions) <(awk '/^## Software versions/ { section = 1 }
                          section && /^```/ { block++; next }
                          section && block == 1 { print }' README.md)
else
  versions
fi
