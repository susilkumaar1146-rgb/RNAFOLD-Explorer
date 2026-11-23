"""Tests for modules/energy_eval.py"""

import pytest

from modules.energy_eval import decompose_loops, evaluate_structure
from modules.energy_parameters import EnergyModel

MODEL = EnergyModel()


def test_decompose_loops_children():
    # ((..((...))..((...)).))   -> outer pair has two children (multiloop)
    pairs = [(0, 22), (1, 20), (4, 10), (5, 9), (13, 19)]
    exterior, children = decompose_loops(pairs, 23)
    assert exterior == [(0, 22)]
    assert children[(0, 22)] == [(1, 20)]
    assert children[(1, 20)] == [(4, 10), (13, 19)]
    assert children[(5, 9)] == []


def test_slide_wild_type_breakdown():
    ev = evaluate_structure("GGGGAAAACCCC", [(0, 11), (1, 10), (2, 9), (3, 8)])
    types = [l.loop_type for l in ev.loops]
    assert types.count("Stack") == 3
    assert types.count("Hairpin") == 1
    assert ev.total == pytest.approx(-4.18)


def test_terminal_au_penalty_in_exterior_and_hairpin():
    # Two A-U pairs (A1-U9, A2-U8) enclosing a 5-nt hairpin loop:
    # the outer A-U ends a helix in the exterior loop, the inner one ends
    # the helix next to the hairpin -> each gets +0.45 kcal/mol.
    seq = "AAGGGGGUU"
    ev = evaluate_structure(seq, [(0, 8), (1, 7)])
    ext = [l for l in ev.loops if l.loop_type == "Exterior"][0]
    hp = [l for l in ev.loops if l.loop_type == "Hairpin"][0]
    assert ext.energy == pytest.approx(0.45)
    assert hp.energy == pytest.approx(MODEL.hairpin_initiation(5) + 0.45)


def test_multiloop_energy_formula():
    seq = "GGGGGAGGGGAAAACCCCAGGGGAAAACCCCACCCCC"
    #      0   4 6   10  14 18  23  27 31  36
    pairs = [(0, 36), (1, 35), (2, 34), (3, 33), (4, 32),
             (6, 17), (7, 16), (8, 15), (9, 14),
             (19, 30), (20, 29), (21, 28), (22, 27)]
    ev = evaluate_structure(seq, pairs)
    ml = [l for l in ev.loops if l.loop_type == "Multiloop"][0]
    # a + 3 helices * c, all G-C ends, unpaired cost b = 0
    assert ml.energy == pytest.approx(MODEL.ml_a + 3 * MODEL.ml_c + ml.unpaired * MODEL.ml_b)
    assert ml.unpaired == 3


def test_invalid_structure_rejected():
    with pytest.raises(ValueError):
        evaluate_structure("GGGGAAAACCCC", [(0, 11), (0, 10)])
    with pytest.raises(ValueError):
        evaluate_structure("AAAAAAAAAAAA", [(0, 11)])     # A-A not a pair


def test_empty_structure_energy_zero():
    assert evaluate_structure("GGGGAAAACCCC", []).total == 0.0
