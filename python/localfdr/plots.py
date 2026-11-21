# -----------------------------------------------------------------------------
# File:    python/localfdr/plots.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Figures 1, S1, S2, S3.
# Purpose: The plotting functions of the figures and the names of their files.
# Inputs:  None.
# Outputs: None.
# Usage:   from localfdr.plots import frequency_barplot, rate_plot
# Source:  R/plots.R. The figures are drawn with matplotlib and show the same
#          marks as those of R (bars, lines, points, tiles, in the same order
#          and colours); the layout around them is matplotlib's.
# -----------------------------------------------------------------------------
import math

import matplotlib
import numpy as np
from matplotlib.backends.backend_pdf import FigureCanvasPdf
from matplotlib.figure import Figure

from . import rbase
from .real_data import DOMAINS

# Names of the six methods in the legend of Figure S1 (drawn in this order).
PROCEDURES = {"NP": "Bayesian - NP", "GP": "Bayesian - P", "SP": "Bayesian - SP",
              "OS": "EB - One Stage", "TS": "EB - Two Stage", "Storey": "Storey"}

# Colours (ColorBrewer Dark2), line types and point shapes of the six
# procedures of Figure S1, in the order of their names.
DARK2 = ("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")
LINE_TYPES = ("solid", "solid", "solid", "dashed", "dashed", "dashed")
POINT_SHAPES = ("*", "s", "^", "o", "^", "+")  # R shapes 8, 15, 17, 1, 2, 3
HOLLOW_SHAPES = (False, False, False, True, True, False)

# Shades of Figure S3: lower ends of the FDR classes, fill and text colours.
SHADE_BREAKS = (0.05, 0.06, 0.07)
SHADE_LABELS = ("below 0.05", "0.05 to 0.06", "0.06 to 0.07", "0.07 and above")
SHADE_FILLS = ("#FCE4D6", "#F4B084", "#C65911", "#833C0B")
SHADE_TEXT = ("black", "black", "black", "white")
MISSING_FILL = "#7F7F7F"

BAR_WIDTH, BAR_SPACE = 1.0, 0.2  # of barplot() in R


def _save(figure, file):
    """Write a figure as PDF, its page content not compressed."""
    with matplotlib.rc_context({"pdf.compression": 0, "pdf.fonttype": 42}):
        FigureCanvasPdf(figure).print_pdf(file)


def _bar_left(n):
    """Left ends of n bars as barplot() places them."""
    return BAR_SPACE + (BAR_WIDTH + BAR_SPACE) * np.arange(n)


def frequency_barplot(frequency, file, size):
    """Bar plot of the number of positions at each count (Figure 1).

    Parameters: frequency (dict from count to frequency,
    `count_frequencies()`), file (path of the PDF to write), size (width and
    height in inches).
    Returns: the plotted frequencies.
    """
    heights = [float(value) for value in frequency.values()]
    figure = Figure(figsize=(size, size))
    axes = figure.add_subplot()
    left = _bar_left(len(heights))
    axes.bar(left, heights, BAR_WIDTH, align="edge", color="#BEBEBE", edgecolor="black",
             linewidth=0.75)
    axes.set_xticks(left + BAR_WIDTH / 2, [str(name) for name in frequency])
    axes.tick_params(axis="x", length=0)
    for side in ("top", "right", "bottom"):
        axes.spines[side].set_visible(False)
    _save(figure, file)
    return frequency


def mixture_barplot(null, alt, file, size):
    """Bar plot of a simulated data set: null counts in gray, non-null counts
    in black, and a dashed line right above the cutoff C = 5 (Figure S2).

    Parameters: null, alt (null and non-null counts), file (path of the PDF
    to write), size (width and height in inches).
    Returns: dict with the plotted frequencies `null` and `alt`.
    """
    K = int(np.max(alt))
    null_frequency = rbase.count_frequencies_from_zero(null, K)
    alt_frequency = rbase.count_frequencies_from_zero(alt, K)
    figure = Figure(figsize=(size, size))
    axes = figure.add_subplot()
    left = _bar_left(K + 1)
    axes.bar(left, null_frequency, BAR_WIDTH, align="edge", color="#BEBEBE", edgecolor="black",
             linewidth=0.75)
    axes.bar(left, alt_frequency, BAR_WIDTH, align="edge", color="#000000", edgecolor="black",
             linewidth=0.75)
    axes.axvline(7.2, linewidth=1.5, linestyle="dashed", color="blue")
    axes.set_xticks(left + BAR_WIDTH / 2, [str(count) for count in range(K + 1)])
    axes.tick_params(axis="x", length=0)
    for side in ("top", "right", "bottom"):
        axes.spines[side].set_visible(False)
    _save(figure, file)
    return {"null": null_frequency, "alt": alt_frequency}


def rate_plot(values, measure, file, size):
    """Line plot of FDR or TPR against N for the six methods (Figure S1).

    Parameters: values (DataFrame with N, method, names of `PROCEDURES`,
    and value), measure ("FDR": with a line at 0.05 and no legend; "TPR":
    with the legend), file (path of the PDF to write), size (width and
    height in inches).
    Returns: the plot (a matplotlib Figure); its lines are the six
    procedures in the order of their names.
    """
    figure = Figure(figsize=(size, size))
    axes = figure.add_subplot()
    procedure = values["method"].map(PROCEDURES)
    for k, name in enumerate(sorted(set(procedure))):
        rows = values[procedure == name].sort_values("N", kind="stable")
        axes.plot(rows["N"].to_numpy(dtype=float), rows["value"].to_numpy(dtype=float),
                  color=DARK2[k], linestyle=LINE_TYPES[k], linewidth=0.75, marker=POINT_SHAPES[k],
                  markersize=4, markerfacecolor="none" if HOLLOW_SHAPES[k] else DARK2[k],
                  label=name)
    axes.set_xlabel("N")
    axes.set_ylabel(measure)
    axes.grid(True, color="#EBEBEB", linewidth=0.5)
    axes.set_axisbelow(True)
    if measure == "FDR":
        axes.axhline(0.05, color="black", linewidth=1.0, linestyle="solid")
    else:
        axes.legend(title="Procedure", loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    figure.tight_layout()
    _save(figure, file)
    return figure


def fdr_shade(value):
    """Class of an FDR value in Figure S3: 0 (below 0.05), 1 (0.05 to 0.06),
    2 (0.06 to 0.07), 3 (0.07 and above); None for a missing value.

    Parameters: value (float). Returns: int, or None.
    """
    if math.isnan(value):
        return None
    return sum(value >= lower for lower in SHADE_BREAKS)


def fdr_table_plot(fdr, methods, file):
    """Shaded table of FDR values by N, lambda, method and theta (Figure S3).

    Cells are shaded by FDR: below 0.05, 0.05 to 0.06, 0.06 to 0.07, 0.07
    and above.

    Parameters: fdr (DataFrame with N, lambda, theta, method and FDR),
    methods (methods in the order of the columns), file (path of the PDF to
    write).
    Returns: the plot (a matplotlib Figure); one rectangle per cell, in the
    order of the rows of `fdr`.
    """
    thetas = sorted(set(fdr["theta"].tolist()))
    columns = [(method, theta) for method in methods for theta in thetas]
    rows = sorted(set(zip(fdr["N"].tolist(), fdr["lambda"].tolist())))
    figure = Figure(figsize=(14, 6.5))
    axes = figure.add_subplot()
    for N, lam, theta, method, value in zip(fdr["N"], fdr["lambda"], fdr["theta"], fdr["method"],
                                            fdr["FDR"].astype(float)):
        x = columns.index((method, theta))
        y = len(rows) - 1 - rows.index((N, lam))  # the first row is at the top
        shade = fdr_shade(value)
        axes.bar(x - 0.5, 1.0, 1.0, bottom=y - 0.5, align="edge", edgecolor="black", linewidth=0.5,
                 color=MISSING_FILL if shade is None else SHADE_FILLS[shade])
        axes.text(x, y, "NA" if shade is None else "%.4f" % value, ha="center", va="center",
                  fontsize=7.4, color="black" if shade is None else SHADE_TEXT[shade])
    axes.set_xlim(-0.5, len(columns) - 0.5)
    axes.set_ylim(-0.5, len(rows) - 0.5)
    axes.xaxis.tick_top()
    axes.set_xticks(range(len(columns)),
                    [f"{method}\ntheta = {rbase.format_number(theta)}" for method, theta in columns])
    axes.set_yticks(range(len(rows)),
                    [f"N = {rbase.format_number(N)}, lambda = {rbase.format_number(lam)}"
                     for N, lam in reversed(rows)])
    axes.tick_params(length=0)
    for side in axes.spines.values():
        side.set_visible(False)
    handles = [matplotlib.patches.Patch(facecolor=fill, edgecolor="black", label=label)
               for fill, label in zip(SHADE_FILLS, SHADE_LABELS)]
    figure.legend(handles=handles, title="FDR", loc="lower center", ncol=4, frameon=False)
    figure.tight_layout(rect=(0, 0.07, 1, 1))
    _save(figure, file)
    return figure


# ----- Names of the figure files -----

# Panels of Figure S1: setting of run0207 and measure.
FIGS1_PANELS = {"setting": (16, 16, 18, 18), "measure": ("FDR", "TPR", "FDR", "TPR")}
# Panels of Figure S2: N and lambda.
FIGS2_PANELS = {"N": (50, 100, 200) * 3, "lambda": (2, 2, 2, 3, 3, 3, 4, 4, 4)}
_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def fig1_file(i):
    """File name of a panel of Figure 1.

    Parameters: i (panel number, from 1; domain `DOMAINS[i - 1]`).
    Returns: file name.
    """
    return "Fig1%s_mutation-counts_%s.pdf" % (_LETTERS[i - 1], DOMAINS[i - 1])


def figS1_file(i):
    """File name of a panel of Figure S1.

    Parameters: i (panel number, from 1; entry of `FIGS1_PANELS`).
    Returns: file name.
    """
    return "FigS1%s_%s_well-separated_setting%d.pdf" % (
        _LETTERS[i - 1], FIGS1_PANELS["measure"][i - 1], FIGS1_PANELS["setting"][i - 1])


def figS2_file(i):
    """File name of a panel of Figure S2.

    Parameters: i (panel number, from 1; entry of `FIGS2_PANELS`).
    Returns: file name.
    """
    return "FigS2%s_simulated-counts_N%d_lambda%d.pdf" % (
        _LETTERS[i - 1], FIGS2_PANELS["N"][i - 1], FIGS2_PANELS["lambda"][i - 1])
