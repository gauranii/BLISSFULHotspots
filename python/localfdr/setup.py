# -----------------------------------------------------------------------------
# File:    python/localfdr/setup.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. All tables and figures.
# Purpose: The settings of the paper, the project paths, the command line
#          arguments of the scripts and a logging helper. Every script in
#          python/scripts/ imports from this file.
# Inputs:  None.
# Outputs: None.
# Usage:   from localfdr.setup import PAPER, PATHS, parse_args
#          Run the scripts from anywhere; the paths are found from this file.
#          Set the environment variable LOCALFDR_OUTPUT to write results/ and
#          figures/ somewhere else than the repository root.
# Source:  R/setup.R.
# -----------------------------------------------------------------------------
import os
import re
import sys
import time

def sourced_file():
    """File of this module, from which the repository root is found (the
    counterpart of `sourced_file()` of R/setup.R).

    Returns: the path of python/localfdr/setup.py.
    """
    return os.path.abspath(__file__)


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(sourced_file()), "..", ".."))

# ----- Settings of the paper, passed to the functions as `settings` -----
PAPER = {
    "alpha": 0.05,          # level of the local fdr and of the FDR
    "replications": 1000,   # simulated data sets per setting
    "sim_draws": 5000,      # stored iterations of the sampler, simulations
    "sim_burn": 1000,       # burn-in iterations, simulations
    "real_draws": 10000,    # stored iterations, protein domains
    "real_burn": 1000,      # burn-in iterations, protein domains
    "real_seed": 1234,      # seed set before each model, protein domains
    # Dirichlet concentration (NP, SP) per model in the protein domain analysis
    "real_gamma": {"NP": 1.5, "SP": 0.02, "P": 0.5, "GP": 0.5},
    "figS2_seed": 1439,     # seed of the simulated data sets of Figure S2
    "figure_size": 7,       # width and height of every figure, inches
}


# ----- Paths -----

def project_paths(root, output=None):
    """Folders and files of the project.

    Parameters: root (repository root), output (folder in which results/
    and figures/ are written; LOCALFDR_OUTPUT, or else the root, by default).
    Returns: dict of paths: `domains`, `archive`, `results`, `figures`.
    """
    output = output or os.environ.get("LOCALFDR_OUTPUT") or root
    return {"domains": os.path.join(root, "data", "domains"),
            "archive": os.path.join(root, "archive"),
            "results": os.path.join(output, "results"),
            "figures": os.path.join(output, "figures")}


PATHS = project_paths(REPO_ROOT)


def output_file(folder, name):
    """Path of an output file, creating its folder when needed.

    Parameters: folder (output folder, for example `paths["results"]`), name
    (file name).
    Returns: the path folder/name.
    """
    os.makedirs(folder, exist_ok=True)
    return os.path.join(folder, name)


# ----- Command line -----

def parse_args(defaults, args=None):
    """Read `--name=value` arguments and bare `--name` flags.

    Every script accepts `--rerun`: recompute by simulation instead of
    reading the archived runs.

    Parameters: defaults (dict of options and their default values; a value
    given on the command line is converted to the type of its default, a
    flag sets a boolean default to True), args (command line arguments).
    Returns: `defaults` with the values given on the command line, and rerun.
    """
    args = sys.argv[1:] if args is None else args
    opts = {"rerun": False, **defaults}
    for arg in args:
        pair = re.fullmatch(r"--([a-z_]+)(=(.*))?", arg)
        if pair is None:
            raise SystemExit("arguments must look like --name=value or --name")
        name, is_flag = pair.group(1), pair.group(2) is None
        if name not in opts:
            raise SystemExit(f"unknown argument --{name}")
        if is_flag != isinstance(opts[name], bool):
            raise SystemExit(f"--{name} " + ("needs a value" if is_flag else "takes no value"))
        opts[name] = True if is_flag else type(opts[name])(pair.group(3))
    return opts


def simulation_options(settings):
    """Default options of a script that can rerun a simulation study.

    Parameters: settings (settings of the paper, `replications`).
    Returns: dict of defaults for `parse_args()`: `runs` (replications to
    run now), `replications` (replications of the study) and `cores`.
    """
    return {"runs": str(settings["replications"]), "replications": settings["replications"],
            "cores": 4}


def parse_index_set(text, count=False):
    """Turn an index text into indices.

    Parameters: text ("3:7" for 3 to 7, "2,5,9", or a single number N),
    count (if True a single number N means 1 to N).
    Returns: list of integers.
    """
    if re.fullmatch(r"[0-9]+:[0-9]+", text):
        first, last = (int(end) for end in text.split(":"))
        step = 1 if last >= first else -1
        return list(range(first, last + step, step))
    index = [int(item) for item in text.split(",")]
    return list(range(1, index[0] + 1)) if count and len(index) == 1 else index


def study_name(name, opts, settings):
    """Name of an output of a simulation study.

    The name carries the number of replications unless it is that of the
    paper.

    Parameters: name (name of the table or figure output, without
    extension), opts (options with `rerun` and `replications`), settings
    (settings of the paper).
    Returns: `name`, or `name_first<D>-replications`.
    """
    if not opts["rerun"] or opts["replications"] == settings["replications"]:
        return name
    return "%s_first%d-replications" % (name, opts["replications"])


def log_message(msg, *values):
    """Print a time-stamped message.

    Parameters: msg (format of the message, with % fields), values (values
    for the format). Returns: None.
    """
    print(time.strftime("%Y-%m-%d %H:%M:%S"), "-", msg % values, flush=True)
