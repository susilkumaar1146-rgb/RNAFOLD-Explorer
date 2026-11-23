"""Tests for modules/sequence_validation.py"""

from modules.sequence_validation import (
    MIN_LENGTH, base_composition, gc_content, parse_fasta, validate_sequence,
)


def test_valid_rna_sequence():
    r = validate_sequence("GGGGAAAACCCC")
    assert r.is_valid
    assert r.sequence == "GGGGAAAACCCC"
    assert r.length == 12
    assert r.errors == []


def test_lowercase_is_converted():
    r = validate_sequence("ggggaaaacccc")
    assert r.is_valid and r.sequence == "GGGGAAAACCCC"


def test_whitespace_and_line_breaks_removed():
    r = validate_sequence("GGGG AAAA\nCC\tCC\r\n")
    assert r.is_valid and r.sequence == "GGGGAAAACCCC"


def test_t_converted_to_u_with_warning():
    r = validate_sequence("GGGGTTTTCCCC")
    assert r.is_valid
    assert r.sequence == "GGGGUUUUCCCC"
    assert r.t_converted == 4
    assert any("converted" in w for w in r.warnings)


def test_invalid_characters_rejected_with_positions():
    r = validate_sequence("GGGXAAAN1CCC")
    assert not r.is_valid
    assert r.sequence == ""
    message = r.errors[0]
    assert "'X' at position 4" in message
    assert "'N' at position 8" in message
    assert "'1' at position 9" in message


def test_invalid_characters_are_not_silently_removed():
    r = validate_sequence("GGGG-AAAACCCC")
    assert not r.is_valid           # gap must cause an error, not be dropped
    assert "alignment gap" in r.errors[0]


def test_empty_sequence():
    for raw in ("", "   \n  ", None):
        r = validate_sequence(raw)
        assert not r.is_valid
        assert "empty" in r.errors[0].lower()


def test_too_short_sequence():
    r = validate_sequence("GAAC")
    assert not r.is_valid
    assert any(str(MIN_LENGTH) in e for e in r.errors)
    assert validate_sequence("GAAAC").is_valid      # exactly the minimum


def test_max_length():
    r = validate_sequence("A" * 50, max_length=40)
    assert not r.is_valid
    assert "exceeds" in r.errors[0]


def test_fasta_header_in_text_box_detected():
    r = validate_sequence(">seq1\nGGGGAAAACCCC")
    assert not r.is_valid
    assert "FASTA" in r.errors[0]


def test_parse_fasta_multiple_records():
    text = """>hairpin test
GGGGAAAA
CCCC
; a comment line
>second
gcgcuuuugcgc
"""
    records = parse_fasta(text)
    assert records == [("hairpin test", "GGGGAAAACCCC"), ("second", "gcgcuuuugcgc")]
    assert validate_sequence(records[1][1]).sequence == "GCGCUUUUGCGC"


def test_parse_fasta_without_header():
    assert parse_fasta("GGGGAAAACCCC\n") == [("unnamed_sequence", "GGGGAAAACCCC")]


def test_gc_content_and_composition():
    assert gc_content("GGGGAAAACCCC") == 100 * 8 / 12
    assert gc_content("") == 0.0
    assert base_composition("GGAUC") == {"A": 1, "U": 1, "G": 2, "C": 1}
