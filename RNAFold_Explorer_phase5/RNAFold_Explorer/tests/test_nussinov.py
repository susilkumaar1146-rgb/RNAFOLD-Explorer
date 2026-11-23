"""Tests for modules/nussinov.py"""

import random

from conftest import enumerate_structures
from modules import nussinov
from modules.dot_bracket import validate_dot_bracket


def test_complementary_hairpin():
    r = nussinov.fold("GGGGAAAACCCC")
    assert r.pair_count == 4
    assert r.dot_bracket == "((((....))))"
    assert r.pairs_1based() == [(1, 12), (2, 11), (3, 10), (4, 9)]
    assert r.energy is None                   # Nussinov never reports energy


def test_wobble_pairs_are_used():
    # Only G-U pairs are possible here.
    r = nussinov.fold("GGGAAAUUU")
    assert r.pair_count == 3            # G1-U9, G2-U8, G3-U7 enclosing AAA
    assert r.dot_bracket == "(((...)))"
    for i, j in r.pairs:
        assert {r.sequence[i], r.sequence[j]} == {"G", "U"}


def test_no_possible_pairs():
    r = nussinov.fold("AAAAAAAAAA")
    assert r.pair_count == 0
    assert r.dot_bracket == ".........."
    assert r.optimal_structure_count == 1


def test_minimum_hairpin_loop_constraint():
    # G and C are separated by only 2 bases -> cannot pair
    assert nussinov.fold("GAACAAAA").pair_count == 0
    # separated by exactly 3 bases -> one pair allowed
    r = nussinov.fold("GAAAC")
    assert r.pair_count == 1 and r.dot_bracket == "(...)"


def test_every_hairpin_has_three_unpaired():
    random.seed(1)
    for _ in range(30):
        seq = "".join(random.choice("AUGC") for _ in range(40))
        r = nussinov.fold(seq)
        for i, j in r.pairs:
            assert j - i - 1 >= 3


def test_no_nucleotide_in_two_pairs_and_valid_dot_bracket():
    random.seed(2)
    for _ in range(30):
        seq = "".join(random.choice("AUGC") for _ in range(60))
        r = nussinov.fold(seq)
        positions = [p for pair in r.pairs for p in pair]
        assert len(positions) == len(set(positions))
        assert validate_dot_bracket(r.dot_bracket, seq) == []
        assert len(r.dot_bracket) == len(seq)


def test_traceback_pair_count_matches_matrix():
    seq = "GGGAAAUCCUUCGGAUUUCCC"
    N = nussinov.fill_matrix(seq)
    pairs = nussinov.traceback(seq, N)
    assert len(pairs) == N[0, len(seq) - 1]


def test_matches_brute_force_on_short_sequences():
    """The DP must find the true maximum, and count co-optimal structures exactly."""
    random.seed(3)
    for _ in range(40):
        seq = "".join(random.choice("AUGC") for _ in range(random.randint(5, 13)))
        all_structures = enumerate_structures(seq)
        best = max(len(s) for s in all_structures)
        n_best = sum(1 for s in all_structures if len(s) == best)
        r = nussinov.fold(seq)
        assert r.pair_count == best, seq
        assert r.optimal_structure_count == n_best, seq


def test_mutant_from_slide_has_ties():
    r = nussinov.fold("GGGGAAAACCCA")
    assert r.pair_count == 3
    assert r.alternative_optima
    assert r.optimal_structure_count > 1
