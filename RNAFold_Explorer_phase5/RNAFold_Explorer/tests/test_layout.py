"""Tests for modules/structure_layout.py (pure geometry, no Plotly needed)."""

import math
import random

import pytest

from modules import thermodynamic
from modules.dot_bracket import dot_bracket_to_pairs
from modules.energy_eval import decompose_loops
from modules.structure_layout import (
    arc_points, chord_points, circle_positions, min_nonbonded_distance, planar_layout,
)

TRNA = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"
CLOVERLEAF = "(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))...."

EXAMPLES = [
    ("GGGGAAAACCCC", "((((....))))"),
    ("GGGGAGGGGAAAACCCCCCCC", "((((.((((....))))))))"),
    ("GGGGAAGGGGAAAACCCCAACCCC", "((((..((((....))))..))))"),
    ("GGGCCCAAGCUCGAAAGAGCAAUGGACGAAAGUCCAUAAGGGCCC",
     "((((((..((((....)))).((((((....))))))..))))))"),
    (TRNA, CLOVERLEAF),
]


def test_arc_endpoints_and_height():
    xs, ys = arc_points(0, 11)
    assert xs[0] == pytest.approx(1) and xs[-1] == pytest.approx(12)   # 1-based positions
    assert ys[0] == pytest.approx(0, abs=1e-9) and ys[-1] == pytest.approx(0, abs=1e-9)
    assert max(ys) == pytest.approx(5.5)                                # radius = (12-1)/2
    _, ys_down = arc_points(0, 11, direction=-1)
    assert min(ys_down) == pytest.approx(-5.5)


def test_circle_positions():
    coords = circle_positions(20)
    assert len(coords) == 20
    for x, y in coords:
        assert math.hypot(x, y) == pytest.approx(1.0)
    # first base is just clockwise of the top, and all positions are distinct
    assert coords[0][0] > 0 and coords[0][1] > 0.9
    assert len({(round(x, 6), round(y, 6)) for x, y in coords}) == 20


def test_chord_points_connect_the_two_bases():
    xs, ys = chord_points((1, 0), (0, 1))
    assert (xs[0], ys[0]) == pytest.approx((1, 0))
    assert (xs[-1], ys[-1]) == pytest.approx((0, 1))


@pytest.mark.parametrize("sequence, structure", EXAMPLES)
def test_planar_pairs_and_backbone_have_unit_length(sequence, structure):
    pairs = dot_bracket_to_pairs(structure)
    coords = planar_layout(len(sequence), pairs)
    assert len(coords) == len(sequence)
    for i, j in pairs:                                       # pair width = 1
        assert math.dist(coords[i], coords[j]) == pytest.approx(1.0)
    exterior, _ = decompose_loops(pairs, len(sequence))
    for i, j in exterior:                                    # backbone inside each domain = 1
        for k in range(i, j):
            assert math.dist(coords[k], coords[k + 1]) == pytest.approx(1.0)


@pytest.mark.parametrize("sequence, structure", EXAMPLES)
def test_planar_layout_has_no_overlaps_for_examples(sequence, structure):
    pairs = dot_bracket_to_pairs(structure)
    coords = planar_layout(len(sequence), pairs)
    assert min_nonbonded_distance(coords, pairs) > 0.9


def test_planar_layout_is_deterministic():
    pairs = dot_bracket_to_pairs(CLOVERLEAF)
    assert planar_layout(len(TRNA), pairs) == planar_layout(len(TRNA), pairs)


def test_planar_layout_predicted_structures():
    random.seed(21)
    for _ in range(5):
        seq = "".join(random.choice("AUGC") for _ in range(100))
        result = thermodynamic.fold(seq)
        coords = planar_layout(len(seq), result.pairs)
        assert all(c is not None for c in coords)
        for i, j in result.pairs:
            assert math.dist(coords[i], coords[j]) == pytest.approx(1.0)


def test_unstructured_sequence_is_a_line():
    coords = planar_layout(5, [])
    assert coords == [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0)]
