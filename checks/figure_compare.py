# -----------------------------------------------------------------------------
# File:    checks/figure_compare.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Compare a figure drawn by R with the same figure drawn by Python
#          (matplotlib) by the marks read from the two PDF files. The two
#          lay their pages out differently, so positions are compared after
#          one straight-line map per axis, fitted over all marks of the
#          figure; what is left over is measured in points of the R page.
#          - bars: number, order and fill of the bars; heights proportional;
#            left ends on one straight line.
#          - lines: number of lines, every vertex.
#          - tiles: number of tiles and the fill of each, tile by tile.
# Usage:   from figure_compare import compare_bars, compare_lines, compare_tiles
#          Each returns (marks in the R file, marks in the Python file,
#          numbers compared, largest difference or None when the marks do not
#          correspond).
# -----------------------------------------------------------------------------
import re

from pdf_marks import content

NUMBER = r"(-?\d+(?:\.\d+)?)"
POSITION_TOLERANCE = 0.02  # points of the R page; R prints two decimals
COLOUR_TOLERANCE = 0.002   # of a colour channel in [0, 1]; R prints three decimals


def r_rectangles(path):
    """Rectangles of an R PDF that are filled and stroked: (x, y, width,
    height, (r, g, b)), in drawing order."""
    fill, found = (0.0, 0.0, 0.0), []
    lines = [line.strip() for line in content(path).split("\n")]
    for i, line in enumerate(lines):
        m = re.fullmatch(f"{NUMBER} {NUMBER} {NUMBER} scn", line)
        if m:
            fill = tuple(float(v) for v in m.groups())
        m = re.fullmatch(f"{NUMBER} {NUMBER} {NUMBER} {NUMBER} re", line)
        if m and i + 1 < len(lines) and lines[i + 1] == "B":
            found.append(tuple(float(v) for v in m.groups()) + (fill,))
    return found


def r_polylines(path, vertices):
    """Stroked paths of an R PDF with this number of vertices, as lists of
    (x, y), in drawing order."""
    found = []
    for text in re.findall(r"((?:[\d.]+ [\d.]+ [ml]\n)+)S", content(path)):
        points = [(float(x), float(y)) for x, y in re.findall(r"([\d.]+) ([\d.]+) [ml]", text)]
        if len(points) == vertices:
            found.append(points)
    return found


def matplotlib_paths(path):
    """Paths of a PDF written by matplotlib without compression: (points,
    closed, paint operator, fill (r, g, b)), in drawing order."""
    data = open(path, "rb").read()
    start = data.index(b"stream\n") + 7
    text = data[start:data.index(b"endstream", start)].decode("latin1")
    tokens = text.split()
    fill, points, closed, found = (0.0, 0.0, 0.0), [], False, []
    for i, token in enumerate(tokens):
        if token == "rg":
            fill = tuple(float(v) for v in tokens[i - 3:i])
        elif token == "g":
            fill = (float(tokens[i - 1]),) * 3
        elif token == "m":
            points, closed = [(float(tokens[i - 2]), float(tokens[i - 1]))], False
        elif token == "l":
            points.append((float(tokens[i - 2]), float(tokens[i - 1])))
        elif token == "h":
            closed = True
        elif token in ("B", "S", "f", "n") and points:
            found.append((points, closed, token, fill))
            points = []
    return found


def matplotlib_rectangles(path):
    """Filled and stroked rectangles of a matplotlib PDF: (x, y, width,
    height, (r, g, b)), in drawing order."""
    found = []
    for points, closed, paint, fill in matplotlib_paths(path):
        if closed and paint == "B" and len(points) == 4:
            xs, ys = [p[0] for p in points], [p[1] for p in points]
            found.append((min(xs), min(ys), max(xs) - min(xs), max(ys) - min(ys), fill))
    return found


def matplotlib_polylines(path, vertices):
    """Stroked open paths of a matplotlib PDF with this number of vertices."""
    return [points for points, closed, paint, _ in matplotlib_paths(path)
            if paint == "S" and not closed and len(points) == vertices]


def fit(u, v, through_origin=False):
    """Largest |a + b u - v| of the least squares line from u to v."""
    n = len(u)
    mu, mv = (0.0, 0.0) if through_origin else (sum(u) / n, sum(v) / n)
    sxx = sum((a - mu) ** 2 for a in u)
    b = sum((a - mu) * (c - mv) for a, c in zip(u, v)) / sxx if sxx else 0.0
    a = mv - b * mu
    return max(abs(a + b * x - y) for x, y in zip(u, v))


def colour_difference(a, b):
    return max(abs(x - y) for x, y in zip(a, b))


def compare_bars(r_pdf, py_pdf):
    """Bars of a bar plot: same number, same fill bar by bar, heights
    proportional and left ends on one line (to POSITION_TOLERANCE)."""
    r, py = r_rectangles(r_pdf), matplotlib_rectangles(py_pdf)
    if len(r) != len(py) or not r:
        return len(r), len(py), 0, None
    if max(colour_difference(a[4], b[4]) for a, b in zip(r, py)) > COLOUR_TOLERANCE:
        return len(r), len(py), 0, None
    heights = fit([b[3] for b in py], [a[3] for a in r], through_origin=True)
    lefts = fit([b[0] for b in py], [a[0] for a in r])
    return len(r), len(py), 3 * len(r), max(heights, lefts)


def compare_lines(r_pdf, py_pdf, vertices):
    """Lines of a line plot: same number, every vertex after one map per axis."""
    r, py = r_polylines(r_pdf, vertices), matplotlib_polylines(py_pdf, vertices)
    if len(r) != len(py) or not r:
        return len(r), len(py), 0, None
    r_points = [p for line in r for p in line]
    py_points = [p for line in py for p in line]
    worst = max(fit([p[k] for p in py_points], [p[k] for p in r_points]) for k in (0, 1))
    return len(r), len(py), 2 * len(r_points), worst


def _tiles(rectangles):
    """The rectangles of the most common size, from the top row down, each
    row from the left."""
    widths = sorted(r[2] for r in rectangles)
    typical = widths[len(widths) // 2]
    tiles = [r for r in rectangles if abs(r[2] - typical) < 0.02 * typical]
    return sorted(tiles, key=lambda r: (-round(r[1], 1), r[0]))


def compare_tiles(r_pdf, py_pdf):
    """Tiles of a shaded table: same number, same fill tile by tile; the
    difference returned is that of the colour channels."""
    r, py = _tiles(r_rectangles(r_pdf)), _tiles(matplotlib_rectangles(py_pdf))
    if len(r) != len(py) or not r:
        return len(r), len(py), 0, None
    return len(r), len(py), 3 * len(r), max(colour_difference(a[4], b[4]) for a, b in zip(r, py))
