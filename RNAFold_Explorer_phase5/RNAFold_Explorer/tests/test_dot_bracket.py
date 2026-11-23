"""Tests for modules/dot_bracket.py"""

import pytest

from modules.dot_bracket import (
    dot_bracket_to_pairs, pair_partner_map, pairs_to_dot_bracket,
    validate_dot_bracket, validate_pairs,
)


def test_pairs_to_dot_bracket_simple_hairpin():
    pairs = [(0, 11), (1, 10), (2, 9), (3, 8)]
    assert pairs_to_dot_bracket(pairs, 12) == "((((....))))"


def test_round_trip():
    structure = "((((...))..((...)).))"
    pairs = dot_bracket_to_pairs(structure)
    assert pairs_to_dot_bracket(pairs, len(structure)) == structure


def test_length_equals_sequence_length():
    db = pairs_to_dot_bracket([(2, 8)], 15)
    assert len(db) == 15 and db.count("(") == db.count(")") == 1


def test_empty_structure():
    assert pairs_to_dot_bracket([], 6) == "......"
    assert dot_bracket_to_pairs("......") == []


def test_unbalanced_brackets_raise():
    with pytest.raises(ValueError):
        dot_bracket_to_pairs("((...)")
    with pytest.raises(ValueError):
        dot_bracket_to_pairs("(...)))")


def test_invalid_character_raises():
    with pytest.raises(ValueError):
        dot_bracket_to_pairs("((..x..))")


def test_duplicate_position_rejected():
    errors = validate_pairs([(0, 8), (0, 10)], 12, check_min_hairpin=False)
    assert any("more than one pair" in e for e in errors)
    with pytest.raises(ValueError):
        pairs_to_dot_bracket([(0, 8), (0, 10)], 12)


def test_pseudoknot_rejected():
    errors = validate_pairs([(0, 6), (3, 10)], 12, check_min_hairpin=False)
    assert any("cross" in e for e in errors)


def test_validate_with_sequence():
    assert validate_dot_bracket("((((....))))", "GGGGAAAACCCC") == []
    # A-A is not a valid pair
    assert validate_dot_bracket("((((....))))", "AGGGAAAACCCA") != []
    # hairpin too short
    assert validate_dot_bracket("(..)", "GAAC") != []
    # wrong length
    assert validate_dot_bracket("((....))", "GGGGAAAACCCC") != []


def test_partner_map():
    assert pair_partner_map([(0, 5)], 7) == [5, -1, -1, -1, -1, 0, -1]
