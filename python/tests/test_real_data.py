# -----------------------------------------------------------------------------
# File:    python/tests/test_real_data.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Tests of real_data.py: the reader on a made-up domain file and on
#          the six domains, the frequencies, and the hotspot counts of short
#          chains.
# Usage:   python -m pytest python/tests/test_real_data.py
# -----------------------------------------------------------------------------
import numpy as np
import pytest

from inputs import SHORT
from localfdr import real_data as rd
from localfdr.setup import PATHS


@pytest.fixture
def folder(tmp_path):
    counts = [0] * 30 + [1] * 14 + [2] * 8 + [3] * 4 + [4, 4, 5, 14, 15, 15, 16, 16, 17, 20]
    lines = [str(len(counts)), "0:0:"] + [
        f"{i}:{c}:" + ",".join(f"GENE{i}|p.A{k}T" for k in range(c))
        for i, c in enumerate(counts, start=1)]
    (tmp_path / "made_up.txt").write_text("\n".join(lines) + "\n")
    return tmp_path


def test_read_domain(folder):
    domain = rd.read_domain("made_up", folder)
    assert domain["N"] == 66 and len(domain["counts"]) == 66
    assert len(domain["with_position_0"]) == 67 and domain["with_position_0"][0] == 0
    assert domain["counts"].dtype == np.int64 and domain["counts"][-1] == 20
    assert rd.DOMAINS[0] == "cd00054" and len(rd.DOMAINS) == 6
    for name, N in zip(rd.DOMAINS, (90, 93, 175, 182, 188, 193)):   # Table 3, column N
        domain = rd.read_domain(name, PATHS["domains"])
        assert domain["N"] == N == len(domain["counts"])


def test_count_frequencies():
    assert rd.count_frequencies([0, 0, 2, 4, 4]) == {0: 2, 1: 0, 2: 1, 3: 0, 4: 2}


def test_domain_rejections(folder):
    x = rd.read_domain("made_up", folder)["counts"]
    np_count = rd.domain_bayes_rejections(x, "NP", SHORT)
    assert np_count == rd.domain_bayes_rejections(x, "NP", SHORT)   # the seed is set inside
    assert 6 <= np_count <= 12
    paths = {**PATHS, "domains": str(folder)}
    table = rd.domain_rejections(["made_up"], ["NP", "GP"], SHORT, paths, cores=2)
    assert list(table.columns) == ["domain", "N", "EB", "NP", "GP"]
    assert table["N"][0] == 66 and table["NP"][0] == np_count
    from localfdr.eb_fdr import eb_rejections
    assert table["EB"][0] == eb_rejections(x) and 5 <= table["EB"][0] <= 12
    with pytest.raises(RuntimeError):
        rd.domain_rejections(["made_up"], ["XX"], SHORT, paths)
