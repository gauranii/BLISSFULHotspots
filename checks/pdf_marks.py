# -----------------------------------------------------------------------------
# File:    checks/pdf_marks.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Read the marks of a PDF file written by R: the rectangles (bars)
#          and the polylines (lines of a line plot).
# Usage:   from pdf_marks import rectangles, polylines, same_marks
# -----------------------------------------------------------------------------
import re
import zlib

TOLERANCE = 0.011  # points; the files print coordinates with two decimals


def content(path):
    """The decompressed page content of a PDF file written by R."""
    data = open(path, 'rb').read()
    text = ''
    for m in re.finditer(rb'/Length (\d+) /Filter /FlateDecode\s*>>\s*stream\r?\n', data):
        text += zlib.decompress(data[m.end():m.end() + int(m.group(1))]).decode('latin1')
    return text


def rectangles(path):
    """Every rectangle (x, y, width, height), in drawing order."""
    return [tuple(map(float, r))
            for r in re.findall(r'([\d.]+) ([\d.]+) ([\d.]+) ([\d.]+) re', content(path))]


def polylines(path, vertices):
    """Every stroked path with this number of vertices, as tuples of numbers
    (x1, y1, x2, y2, ...), in drawing order."""
    lines = []
    for path_text in re.findall(r'((?:[\d.]+ [\d.]+ [ml]\n)+)S', content(path)):
        points = re.findall(r'([\d.]+) ([\d.]+) [ml]', path_text)
        if len(points) == vertices:
            lines.append(tuple(float(v) for point in points for v in point))
    return lines


def same_marks(a, b):
    """(number of marks in a, in b, largest coordinate difference or None)."""
    if len(a) != len(b) or any(len(x) != len(y) for x, y in zip(a, b)):
        return len(a), len(b), None
    return len(a), len(b), max((abs(p - q) for x, y in zip(a, b) for p, q in zip(x, y)), default=0.0)
