# -----------------------------------------------------------------------------
# File:    python/tests/test_plots.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Tests of plots.py: every figure is written, shows the marks it should
#          and returns what it plotted; the names of the figure files.
# Usage:   python -m pytest python/tests/test_plots.py
# -----------------------------------------------------------------------------
import math

import numpy as np
import pandas as pd

from localfdr import plots


def bars(figure_or_axes):
    return [patch for patch in figure_or_axes.patches]


def test_bar_positions_and_save(tmp_path):
    np.testing.assert_allclose(plots._bar_left(3), [0.2, 1.4, 2.6])
    from matplotlib.figure import Figure
    plots._save(Figure(), tmp_path / "empty.pdf")
    assert (tmp_path / "empty.pdf").read_bytes().startswith(b"%PDF")


def test_frequency_barplot(tmp_path):
    frequency = {0: 5, 1: 3, 2: 0, 3: 1}
    file = tmp_path / "bars.pdf"
    assert plots.frequency_barplot(frequency, file, 7) == frequency
    content = file.read_bytes()
    assert content.startswith(b"%PDF") and b"FlateDecode" not in content


def test_mixture_barplot(tmp_path):
    plotted = plots.mixture_barplot(np.array([0, 0, 1, 2, 2, 9]), np.array([6, 7, 7]),
                                    tmp_path / "mixture.pdf", 7)
    np.testing.assert_array_equal(plotted["null"], [2, 1, 2, 0, 0, 0, 0, 0])   # 9 is beyond max(alt)
    np.testing.assert_array_equal(plotted["alt"], [0, 0, 0, 0, 0, 0, 1, 2])
    assert (tmp_path / "mixture.pdf").stat().st_size > 1000


def test_rate_plot(tmp_path):
    values = pd.DataFrame({"N": [50, 100, 200] * 6,
                           "method": np.repeat(["NP", "GP", "SP", "TS", "OS", "Storey"], 3),
                           "value": np.linspace(0.01, 0.18, 18)})
    figure = plots.rate_plot(values, "FDR", tmp_path / "fdr.pdf", 7)
    axes = figure.axes[0]
    assert len(axes.lines) == 7 and axes.get_legend() is None   # six procedures and the 0.05 line
    assert [line.get_label() for line in axes.lines[:6]] == sorted(plots.PROCEDURES.values())
    np.testing.assert_array_equal(axes.lines[3].get_ydata(), values["value"][12:15])  # EB - One Stage
    assert axes.lines[3].get_linestyle() == "--" and axes.lines[0].get_linestyle() == "-"
    assert axes.lines[6].get_ydata()[0] == 0.05 and axes.get_ylabel() == "FDR"
    figure = plots.rate_plot(values, "TPR", tmp_path / "tpr.pdf", 7)
    assert len(figure.axes[0].lines) == 6 and figure.axes[0].get_legend() is not None
    assert len(plots.DARK2) == len(plots.LINE_TYPES) == len(plots.POINT_SHAPES) == 6


def test_fdr_shade_and_table_plot(tmp_path):
    assert [plots.fdr_shade(v) for v in (0.0, 0.0499, 0.05, 0.0599, 0.06, 0.07, 0.3)] == \
        [0, 0, 1, 1, 2, 3, 3]
    assert plots.fdr_shade(math.nan) is None
    fdr = pd.DataFrame({"N": [50, 50, 100, 100], "lambda": [1, 1, 1, 1], "theta": [0.15, 0.3] * 2,
                        "method": ["NP"] * 4, "FDR": [0.01, 0.055, 0.065, np.nan]})
    figure = plots.fdr_table_plot(fdr, ["NP"], tmp_path / "table.pdf")
    axes = figure.axes[0]
    assert len(axes.patches) == 4
    from matplotlib.colors import to_hex
    fills = [to_hex(patch.get_facecolor()).upper() for patch in axes.patches]
    assert fills == ["#FCE4D6", "#F4B084", "#C65911", plots.MISSING_FILL]
    assert [text.get_text() for text in axes.texts] == ["0.0100", "0.0550", "0.0650", "NA"]
    assert axes.patches[0].get_y() == 0.5 and axes.patches[2].get_y() == -0.5   # N = 50 on top
    assert [label.get_text() for label in axes.get_yticklabels()] == [
        "N = 100, lambda = 1", "N = 50, lambda = 1"]
    assert axes.get_xticklabels()[1].get_text() == "NP\ntheta = 0.3"


def test_file_names():
    assert plots.fig1_file(1) == "Fig1A_mutation-counts_cd00054.pdf"
    assert plots.fig1_file(6) == "Fig1F_mutation-counts_cd00063.pdf"
    assert plots.figS1_file(2) == "FigS1B_TPR_well-separated_setting16.pdf"
    assert plots.figS1_file(3) == "FigS1C_FDR_well-separated_setting18.pdf"
    assert plots.figS2_file(1) == "FigS2A_simulated-counts_N50_lambda2.pdf"
    assert plots.figS2_file(9) == "FigS2I_simulated-counts_N200_lambda4.pdf"
    assert len(plots.FIGS2_PANELS["N"]) == 9 and len(plots.FIGS1_PANELS["setting"]) == 4
