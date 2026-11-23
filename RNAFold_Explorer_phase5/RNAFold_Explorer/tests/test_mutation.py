"""Tests for modules/mutation_analysis.py"""

import random

import pytest

from modules import thermodynamic
from modules.mutation_analysis import (
    SCAN_MAX_LENGTH, analyse_mutation, compare_pairs, create_mutant, interpret,
    parse_mutation, scan_mutations,
)

WT = "GGGGAAAACCCC"          # ((((....))))  pairs 1-12, 2-11, 3-10, 4-9


# --------------------------- creating mutants --------------------------- #
def test_mutant_sequence_generated_automatically():
    mutant, m = create_mutant(WT, 12, "A")
    assert mutant == "GGGGAAAACCCA"
    assert m.label == "C12A" and m.original == "C" and m.index == 11


def test_mutation_at_first_position():
    mutant, m = create_mutant(WT, 1, "A")
    assert mutant == "AGGGAAAACCCC" and m.label == "G1A"


def test_mutation_at_last_position():
    mutant, _ = create_mutant(WT, len(WT), "G")
    assert mutant == "GGGGAAAACCCG"


def test_t_is_accepted_as_u():
    mutant, m = create_mutant(WT, 5, "t")
    assert mutant[4] == "U" and m.replacement == "U"


@pytest.mark.parametrize("position, base", [(0, "A"), (13, "A"), (1, "G"), (3, "X"), (3, "AU")])
def test_invalid_mutations_rejected(position, base):
    with pytest.raises(ValueError):
        create_mutant(WT, position, base)


def test_parse_mutation():
    m = parse_mutation("c12a", WT)
    assert (m.original, m.position, m.replacement) == ("C", 12, "A")
    with pytest.raises(ValueError):
        parse_mutation("G12A", WT)        # position 12 is C, not G
    with pytest.raises(ValueError):
        parse_mutation("12A")


# ---------------------------- pair comparison --------------------------- #
def test_compare_pairs_sets():
    wt = [(0, 11), (1, 10), (2, 9)]
    mut = [(1, 10), (2, 9), (3, 8)]
    c = compare_pairs(wt, mut, 12)
    assert c.retained == [(1, 10), (2, 9)]
    assert c.lost == [(0, 11)]
    assert c.new == [(3, 8)]
    assert c.changed_positions == [0, 3, 8, 11]
    assert not c.identical


def test_retained_lost_new_partition_both_structures():
    random.seed(5)
    for _ in range(20):
        seq = "".join(random.choice("AUGC") for _ in range(40))
        pos = random.randint(1, 40)
        base = random.choice([b for b in "AUGC" if b != seq[pos - 1]])
        a = analyse_mutation(seq, pos, base)
        c = a.comparison
        assert sorted(c.retained + c.lost) == a.wild_type.pairs
        assert sorted(c.retained + c.new) == a.mutant.pairs


# --------------------------- full analysis ------------------------------ #
def test_slide_mutation_paired_nucleotide():
    a = analyse_mutation(WT, 12, "A", "Thermodynamic")
    assert a.mutant.sequence == "GGGGAAAACCCA"
    assert a.wild_type.pair_count == 4 and a.mutant.pair_count == 3
    assert a.comparison.lost == [(0, 11)]               # G1-C12 (0-based)
    assert a.comparison.new == []
    assert len(a.comparison.retained) == 3
    assert a.wt_energy == pytest.approx(-4.18)
    assert a.mut_energy == pytest.approx(-0.92)
    assert a.ddg == pytest.approx(3.26)
    assert "paired with G1" in a.wt_site
    assert a.mut_site.startswith("unpaired")


def test_ddg_is_mutant_minus_wild_type():
    random.seed(6)
    for _ in range(10):
        seq = "".join(random.choice("AUGC") for _ in range(35))
        a = analyse_mutation(seq, 10, "G" if seq[9] != "G" else "C")
        expected = thermodynamic.fold(a.mutant.sequence).energy - thermodynamic.fold(seq).energy
        assert a.ddg == pytest.approx(round(expected, 2))


def test_mutation_creating_a_new_pair():
    a = analyse_mutation("GGGGAAAACCCA", 12, "C")       # reverse mutation A12C
    assert a.comparison.new == [(0, 11)]
    assert a.ddg == pytest.approx(-3.26)


def test_unpaired_nucleotide_without_structural_change():
    a = analyse_mutation(WT, 6, "U")                    # A6 in the hairpin loop
    assert a.wt_site.startswith("unpaired (hairpin loop)")
    assert a.comparison.identical
    assert a.comparison.changed_positions == []
    assert a.ddg == 0.0
    text = " ".join(interpret(a))
    assert "same base pairs" in text and "no change" in text


def test_wobble_compatible_mutation_keeps_pair_possible():
    # C10 pairs with G3; C -> U leaves a G-U wobble pair possible
    a = analyse_mutation(WT, 10, "U")
    assert "can still pair" in interpret(a)[0]


def test_nussinov_mutation_analysis():
    a = analyse_mutation(WT, 12, "A", "Nussinov")
    assert a.wild_type.pair_count == 4 and a.mutant.pair_count == 3
    assert a.wild_type.energy is None                   # Nussinov reports no MFE
    assert "Nussinov does not use energies" in a.energy_note
    assert any("ties are common" in line for line in interpret(a))


def test_wild_type_result_reused():
    wt = thermodynamic.fold(WT)
    a = analyse_mutation(WT, 12, "A", wild_type=wt)
    assert a.wild_type is wt


def test_interpretation_makes_no_functional_claim():
    a = analyse_mutation(WT, 12, "A")
    text = " ".join(interpret(a)).lower()
    assert "does not by itself show a loss or gain of biological function" in text
    assert "predicted to be less stable under this model" in text


def test_summary_rows():
    rows = analyse_mutation(WT, 12, "A").summary_rows()
    values = {r["Quantity"]: r["Value"] for r in rows}
    assert values["Mutation"] == "C12A"
    assert values["ΔΔG (mutant − wild type)"] == "+3.26"


# ------------------------------- scan ----------------------------------- #
def test_scan_covers_every_single_mutant():
    rows = scan_mutations(WT)
    assert len(rows) == 3 * len(WT)
    labels = {r["Mutation"] for r in rows}
    assert "C12A" in labels and "A6U" in labels
    c12a = [r for r in rows if r["Mutation"] == "C12A"][0]
    assert c12a["ΔΔG (kcal/mol)"] == pytest.approx(3.26)
    assert c12a["WT site paired"] is True


def test_scan_progress_and_limit():
    seen = []
    scan_mutations("GGGAAACCC", progress=seen.append)
    assert seen[-1] == pytest.approx(1.0)
    with pytest.raises(ValueError):
        scan_mutations("A" * (SCAN_MAX_LENGTH + 1))
