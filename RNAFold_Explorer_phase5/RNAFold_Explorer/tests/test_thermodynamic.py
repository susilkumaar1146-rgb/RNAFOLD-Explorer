"""Tests for modules/thermodynamic.py and modules/energy_parameters.py"""

import itertools
import math
import random

import pytest

from conftest import enumerate_structures
from modules import thermodynamic
from modules.dot_bracket import validate_dot_bracket
from modules.energy_eval import evaluate_structure
from modules.energy_parameters import EnergyModel
from modules.pairing import ALLOWED_PAIRS

MODEL = EnergyModel()


# --------------------------- parameter table --------------------------- #
def test_stacking_table_covers_all_pair_combinations():
    """Every combination of two allowed pairs (6 x 6 = 36) needs a stack value."""
    for (x, a), (y, b) in itertools.product(ALLOWED_PAIRS, repeat=2):
        assert f"{x}{y}/{a}{b}" in MODEL.stack_table


def test_published_values_loaded():
    assert MODEL.stack_table["GC/CG"] == -3.42
    assert MODEL.stack_table["GG/CC"] == -3.26
    assert MODEL.stack_table["CC/GG"] == -3.26           # symmetric equivalent
    assert MODEL.hairpin_initiation(4) == 5.6
    assert MODEL.terminal_penalty_value == 0.45


def test_hairpin_extrapolation_matches_table_trend():
    # The 1.75 RT ln(n/9) formula reproduces the tabulated value at n = 30.
    g30 = 6.4 + 1.75 * MODEL.RT * math.log(30 / 9)
    assert abs(g30 - MODEL.hairpin_initiation(30)) < 0.05
    assert MODEL.hairpin_initiation(40) > MODEL.hairpin_initiation(30)


# ------------------------- worked example (slide 5) --------------------- #
def test_slide_wild_type():
    r = thermodynamic.fold("GGGGAAAACCCC")
    assert r.dot_bracket == "((((....))))"
    # 3 x GG/CC (-3.26) + hairpin(4) 5.6 = -4.18
    assert r.energy == pytest.approx(-4.18, abs=1e-9)


def test_slide_mutant_c12a():
    r = thermodynamic.fold("GGGGAAAACCCA")
    assert r.dot_bracket == ".(((....)))."
    assert r.energy == pytest.approx(-0.92, abs=1e-9)
    ddg = r.energy - thermodynamic.fold("GGGGAAAACCCC").energy
    assert ddg == pytest.approx(3.26, abs=1e-9)


# ------------------------------ loop types ------------------------------ #
def loop_types(sequence, result):
    return [loop.loop_type for loop in evaluate_structure(sequence, result.pairs).loops]


def test_hairpin_formation():
    r = thermodynamic.fold("GGGAAACCC")
    assert r.dot_bracket == "(((...)))"
    assert r.energy == pytest.approx(2 * -3.26 + 5.4)
    assert "Hairpin" in loop_types("GGGAAACCC", r)


def test_stacking_sum_is_correct():
    seq = "GCGCAAAAGCGC"
    ev = evaluate_structure(seq, [(0, 11), (1, 10), (2, 9), (3, 8)])
    stacks = [l.energy for l in ev.loops if l.loop_type == "Stack"]
    assert stacks == [MODEL.stack_table["GC/CG"], MODEL.stack_table["CG/GC"], MODEL.stack_table["GC/CG"]]


def test_single_nucleotide_bulge_predicted():
    seq = "GGGGAGGGGAAAACCCCCCCC"
    r = thermodynamic.fold(seq)
    assert "Bulge" in loop_types(seq, r)
    bulge = [l for l in evaluate_structure(seq, r.pairs).loops if l.loop_type == "Bulge"][0]
    # 1-nt bulge = initiation(1) + stack of the two closing pairs (GG/CC)
    assert bulge.energy == pytest.approx(3.8 + -3.26)


def test_larger_bulge_energy_rule():
    # outer pair G1-C10, inner pair G4-C9: 2 unpaired (A2, A3) on the 5' side only
    assert MODEL.interior("GAAGAAAACC", 0, 9, 3, 8) == pytest.approx(2.8)   # G-C ends
    # same bulge closed by an A-U outer pair: + one terminal AU penalty
    assert MODEL.interior("AAAGAAAACU", 0, 9, 3, 8) == pytest.approx(2.8 + 0.45)
    assert MODEL.interior_loop_type(0, 9, 3, 8) == "Bulge"


def test_internal_loop_predicted():
    seq = "GGGGAAGGGGAAAACCCCAACCCC"
    r = thermodynamic.fold(seq)
    assert "Internal loop" in loop_types(seq, r)


def test_internal_loop_energy_rule():
    # 2x2 internal loop between G-C pairs: init(4) = 1.1, symmetric, no AU/GU
    assert MODEL.interior("GAAGAAAACAAC", 0, 11, 3, 8) == pytest.approx(1.1)
    # 1x3 loop closed by A-U and G-C: init(4) + 0.6*2 + 0.7
    assert MODEL.interior("AAGAAAACAAAU", 0, 11, 2, 7) == pytest.approx(1.1 + 1.2 + 0.7)


def test_multiloop_predicted():
    # designed: outer 6-bp stem enclosing two hairpins
    seq = "GGGCCCAAGCUCGAAAGAGCAAUGGACGAAAGUCCAUAAGGGCCC"
    r = thermodynamic.fold(seq)
    assert "Multiloop" in loop_types(seq, r)
    assert validate_dot_bracket(r.dot_bracket, seq) == []


def test_unstructured_sequence():
    r = thermodynamic.fold("AAAAAAAAAA")
    assert r.pairs == [] and r.energy == 0.0
    assert r.dot_bracket == ".........."


# --------------------- consistency and optimality ----------------------- #
def test_reported_energy_equals_structure_energy():
    random.seed(10)
    for _ in range(20):
        seq = "".join(random.choice("AUGC") for _ in range(70))
        r = thermodynamic.fold(seq)
        ev = evaluate_structure(seq, r.pairs)
        assert r.energy == pytest.approx(round(ev.total, 2), abs=1e-9)
        assert validate_dot_bracket(r.dot_bracket, seq) == []


def test_dp_finds_true_minimum_on_short_sequences():
    """Brute force: evaluate every possible structure and compare."""
    random.seed(11)
    for _ in range(40):
        seq = "".join(random.choice("AUGC") for _ in range(random.randint(5, 14)))
        best = min(evaluate_structure(seq, s).total for s in enumerate_structures(seq))
        assert thermodynamic.fold(seq).energy == pytest.approx(round(best, 2), abs=1e-9), seq


def test_dp_minimum_gc_rich_sequences():
    """GC-rich sequences have many structures with multiloops and bulges."""
    random.seed(12)
    for _ in range(15):
        seq = "".join(random.choice("GGCCAU") for _ in range(16))
        best = min(evaluate_structure(seq, s).total for s in enumerate_structures(seq))
        assert thermodynamic.fold(seq).energy == pytest.approx(round(best, 2), abs=1e-9), seq


def test_length_limit():
    with pytest.raises(ValueError):
        thermodynamic.fold("A" * (thermodynamic.MAX_LENGTH + 1))
