"""
structure_layout.py
===================
Computes WHERE to draw each nucleotide and base pair. Pure Python (no
plotting library), so the geometry can be unit-tested on its own;
visualization.py turns these coordinates into Plotly figures.

Three layouts
-------------
1. Arc diagram  - bases on a horizontal line (x = position), each pair drawn
                  as a semicircular arc. Nested pairs give nested arcs, so the
                  stems are easy to see and the picture is always exact.
2. Circle plot  - bases evenly spaced on a circle, each pair drawn as a curved
                  chord. Also exact; nested structures never cross.
3. Planar layout- a SCHEMATIC 2-D drawing, like textbook RNA figures:
                  * every stem is a straight ladder (pairs 1 unit apart,
                    pair width 1 unit);
                  * every loop is a regular polygon whose corners are the
                    loop's unpaired bases and the bases of its closing pairs;
                  * branches leave each loop perpendicular to the loop.
                  This is deterministic (same input -> same picture) and is
                  NOT a three-dimensional or experimentally determined shape.

Planar layout algorithm (recursive "loop polygon" method)
---------------------------------------------------------
place_stem(i, j):  put pair (i, j); while the next pair (i+1, j-1) is stacked
                   on it, put it one unit further along the stem direction.
place_loop(a, b):  the last pair (a, b) of a stem closes a loop with m corners
                   (a, the loop's unpaired bases, both bases of each inner pair,
                   b). All edges are 1 unit long, so the corners lie on a circle
                   of radius R = 1 / (2 sin(pi/m)). Each inner pair starts a new
                   stem pointing away from the circle's centre -> recursion.
Exterior loop:     top-level stems are laid out separately, then placed
                   side by side along the x-axis so they do not overlap.
Large or highly branched structures can still show some overlap - an
accepted limitation of simple schematic layouts.
"""

import math
from typing import Dict, List, Optional, Tuple

from modules.energy_eval import decompose_loops

Pair = Tuple[int, int]
Point = Tuple[float, float]


# ------------------------------------------------------------------ #
# 1. Arc diagram
# ------------------------------------------------------------------ #
def arc_points(i: int, j: int, direction: int = 1, n_points: int = 40) -> Tuple[List[float], List[float]]:
    """
    Semicircle joining 1-based x positions i+1 and j+1 (inputs are 0-based).
    direction = +1 draws above the axis, -1 below (used for mutant comparison).
    """
    x1, x2 = i + 1, j + 1
    centre, radius = (x1 + x2) / 2, (x2 - x1) / 2
    xs, ys = [], []
    for k in range(n_points + 1):
        t = math.pi - math.pi * k / n_points          # from left end to right end
        xs.append(centre + radius * math.cos(t))
        ys.append(direction * radius * math.sin(t))
    return xs, ys


# ------------------------------------------------------------------ #
# 2. Circle plot
# ------------------------------------------------------------------ #
def circle_positions(n: int, radius: float = 1.0) -> List[Point]:
    """
    Bases clockwise from the top of the circle. One empty slot is left between
    the 3' end and the 5' end so the two ends are visibly separate.
    """
    coords = []
    for k in range(n):
        angle = math.pi / 2 - 2 * math.pi * (k + 0.5) / (n + 1)
        coords.append((radius * math.cos(angle), radius * math.sin(angle)))
    return coords


def chord_points(p1: Point, p2: Point, bend: float = 0.35, n_points: int = 30) -> Tuple[List[float], List[float]]:
    """
    Quadratic Bezier curve from p1 to p2 whose control point is the chord's
    midpoint pulled towards the centre (bend = 0 gives a straight line).
    """
    cx, cy = (p1[0] + p2[0]) / 2 * bend, (p1[1] + p2[1]) / 2 * bend
    xs, ys = [], []
    for k in range(n_points + 1):
        t = k / n_points
        a, b, c = (1 - t) ** 2, 2 * (1 - t) * t, t ** 2
        xs.append(a * p1[0] + b * cx + c * p2[0])
        ys.append(a * p1[1] + b * cy + c * p2[1])
    return xs, ys


# ------------------------------------------------------------------ #
# 3. Planar (schematic 2-D) layout
# ------------------------------------------------------------------ #
def _unit(vx: float, vy: float) -> Point:
    length = math.hypot(vx, vy)
    return vx / length, vy / length


def planar_layout(n: int, pairs: List[Pair]) -> List[Point]:
    """Return one (x, y) coordinate per nucleotide (0-based index)."""
    pairs = sorted(pairs)
    exterior, children = decompose_loops(pairs, n)
    coords: List[Optional[Point]] = [None] * n

    def place_stem(i: int, j: int, pi: Point, pj: Point, direction: Point) -> None:
        coords[i], coords[j] = pi, pj
        while children[(i, j)] == [(i + 1, j - 1)]:          # stacked pair
            pi = (pi[0] + direction[0], pi[1] + direction[1])
            pj = (pj[0] + direction[0], pj[1] + direction[1])
            i, j = i + 1, j - 1
            coords[i], coords[j] = pi, pj
        place_loop(i, j, direction)

    def place_loop(a: int, b: int, direction: Point) -> None:
        # Corners of the loop polygon, in 5'->3' order.
        corners: List[int] = [a]
        inner = {p: q for p, q in children[(a, b)]}
        k = a + 1
        while k < b:
            corners.append(k)
            if k in inner:
                corners.append(inner[k])
                k = inner[k] + 1
            else:
                k += 1
        corners.append(b)
        m = len(corners)

        radius = 1 / (2 * math.sin(math.pi / m))
        apothem = radius * math.cos(math.pi / m)          # centre-to-edge distance
        pa, pb = coords[a], coords[b]
        mid = ((pa[0] + pb[0]) / 2, (pa[1] + pb[1]) / 2)
        centre = (mid[0] + direction[0] * apothem, mid[1] + direction[1] * apothem)

        # Walk around the circle in the direction that keeps a and b fixed.
        ux, uy = _unit(pb[0] - pa[0], pb[1] - pa[1])
        orientation = 1 if (ux * direction[1] - uy * direction[0]) > 0 else -1
        step = -orientation * 2 * math.pi / m
        start = math.atan2(pa[1] - centre[1], pa[0] - centre[0])
        for idx in range(1, m - 1):
            angle = start + idx * step
            coords[corners[idx]] = (centre[0] + radius * math.cos(angle),
                                    centre[1] + radius * math.sin(angle))

        for p, q in children[(a, b)]:                     # branches point outwards
            pp, pq = coords[p], coords[q]
            out = _unit((pp[0] + pq[0]) / 2 - centre[0], (pp[1] + pq[1]) / 2 - centre[1])
            place_stem(p, q, pp, pq, out)

    # Exterior loop: lay out each top-level stem, then shift it into place.
    cursor = 0.0
    right_edge = -math.inf
    k = 0
    top_level = {i: j for i, j in exterior}
    while k < n:
        if k in top_level:
            j = top_level[k]
            place_stem(k, j, (0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
            xs = [coords[t][0] for t in range(k, j + 1)]
            offset = max(cursor, right_edge + 1.5 - min(xs))
            for t in range(k, j + 1):
                coords[t] = (coords[t][0] + offset, coords[t][1])
            right_edge = offset + max(xs)
            cursor = offset + 2.0
            k = j + 1
        else:
            coords[k] = (cursor, 0.0)
            right_edge = max(right_edge, cursor)
            cursor += 1.0
            k += 1
    return [c if c is not None else (0.0, 0.0) for c in coords]


def min_nonbonded_distance(coords: List[Point], pairs: List[Pair]) -> float:
    """
    Smallest distance between two nucleotides that are neither backbone
    neighbours nor base-paired. Values well below 1 indicate overlapping
    parts of the drawing (used by the tests to check layout quality).
    """
    bonded = {(i, j) for i, j in pairs}
    best = math.inf
    n = len(coords)
    for a in range(n):
        for b in range(a + 2, n):
            if (a, b) in bonded:
                continue
            d = math.dist(coords[a], coords[b])
            best = min(best, d)
    return best


def layout_bounds(coords: List[Point]) -> Dict[str, float]:
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    return {"xmin": min(xs), "xmax": max(xs), "ymin": min(ys), "ymax": max(ys)}
