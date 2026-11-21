# -----------------------------------------------------------------------------
# File:    python/tests/test_setup.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Tests of setup.py: the settings of the paper, the paths, the command
#          line arguments and the names of the outputs.
# Usage:   python -m pytest python/tests/test_setup.py
# -----------------------------------------------------------------------------
import os
import re

import pytest

from localfdr import setup


def test_settings_and_paths(tmp_path, monkeypatch):
    assert setup.PAPER["alpha"] == 0.05 and setup.PAPER["real_seed"] == 1234
    assert setup.PAPER["real_gamma"] == {"NP": 1.5, "SP": 0.02, "P": 0.5, "GP": 0.5}
    assert setup.sourced_file().endswith(os.path.join("python", "localfdr", "setup.py"))
    assert os.path.isdir(os.path.join(setup.REPO_ROOT, "R"))
    paths = setup.project_paths("root", "out")
    assert paths == {"domains": os.path.join("root", "data", "domains"),
                     "archive": os.path.join("root", "archive"),
                     "results": os.path.join("out", "results"),
                     "figures": os.path.join("out", "figures")}
    monkeypatch.setenv("LOCALFDR_OUTPUT", "elsewhere")
    assert setup.project_paths("root")["results"] == os.path.join("elsewhere", "results")
    monkeypatch.delenv("LOCALFDR_OUTPUT")
    assert setup.project_paths("root")["figures"] == os.path.join("root", "figures")
    file = setup.output_file(str(tmp_path / "new" / "folder"), "a.csv")
    assert os.path.isdir(os.path.dirname(file)) and file.endswith("a.csv")


def test_parse_args():
    defaults = {"cores": 4, "domains": "all", "poisson": False}
    assert setup.parse_args(defaults, []) == {"rerun": False, **defaults}
    opts = setup.parse_args(defaults, ["--cores=2", "--poisson", "--domains=cd00083,cd00054",
                                       "--rerun"])
    assert opts == {"rerun": True, "cores": 2, "domains": "cd00083,cd00054", "poisson": True}
    for wrong in (["cores=2"], ["--unknown=1"], ["--cores"], ["--poisson=yes"]):
        with pytest.raises(SystemExit):
            setup.parse_args(defaults, wrong)
    assert setup.simulation_options(setup.PAPER) == {"runs": "1000", "replications": 1000,
                                                     "cores": 4}


def test_parse_index_set_and_study_name():
    assert setup.parse_index_set("3:7") == [3, 4, 5, 6, 7]
    assert setup.parse_index_set("2,5,9") == [2, 5, 9]
    assert setup.parse_index_set("4") == [4] and setup.parse_index_set("4", count=True) == [1, 2, 3, 4]
    assert setup.parse_index_set("3:1") == [3, 2, 1]
    paper = setup.PAPER
    assert setup.study_name("Table1", {"rerun": False, "replications": 5}, paper) == "Table1"
    assert setup.study_name("Table1", {"rerun": True, "replications": 1000}, paper) == "Table1"
    assert setup.study_name("Table1", {"rerun": True, "replications": 5}, paper) == \
        "Table1_first5-replications"


def test_log_message(capsys):
    setup.log_message("wrote %d files to %s", 3, "results")
    assert re.fullmatch(r"\d{4}-\d\d-\d\d \d\d:\d\d:\d\d - wrote 3 files to results\n",
                        capsys.readouterr().out)
