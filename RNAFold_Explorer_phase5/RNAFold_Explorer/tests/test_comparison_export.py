"""Tests for modules/comparison.py and modules/results_export.py"""

import io
import random
import zipfile

import pytest

from modules import nussinov, thermodynamic
from modules.comparison import compare_predictions
from modules.mutation_analysis import analyse_mutation, scan_mutations
from modules.results_export import (
    aligned_block, bundle_zip, dot_bracket_text, fasta, mutation_csv, pairs_csv,
    safe_name, scan_csv, summary_report,
)

TRNA = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"


# ------------------------------ comparison ------------------------------ #
def test_identical_predictions():
    seq = "GGGGAAAACCCC"
    c = compare_predictions(nussinov.fold(seq), thermodynamic.fold(seq))
    assert c.identical and c.position_agreement == 100.0
    assert c.nussinov_energy == pytest.approx(-4.18)
    assert "exactly the same" in c.explanation()[0]


def test_trna_predictions_differ():
    c = compare_predictions(nussinov.fold(TRNA), thermodynamic.fold(TRNA))
    assert not c.identical
    assert c.nussinov.pair_count >= c.thermodynamic.pair_count
    assert c.nussinov_energy > c.thermodynamic.energy
    assert sorted(c.shared + c.only_nussinov) == c.nussinov.pairs
    assert sorted(c.shared + c.only_thermodynamic) == c.thermodynamic.pairs
    rows = {r["Parameter"]: r for r in c.table_rows()}
    assert rows["Estimated free energy"]["Nussinov"] == "Not calculated by Nussinov"
    assert rows["Stacking energy"]["Nussinov"] == "Not considered"


def test_comparison_invariants_random():
    """Nussinov never has fewer pairs; thermodynamic never has higher energy (short loops)."""
    random.seed(30)
    for _ in range(15):
        seq = "".join(random.choice("AUGC") for _ in range(50))
        c = compare_predictions(nussinov.fold(seq), thermodynamic.fold(seq))
        assert c.nussinov.pair_count >= c.thermodynamic.pair_count
        assert c.thermodynamic.energy <= c.nussinov_energy + 1e-9


def test_different_sequences_rejected():
    with pytest.raises(ValueError):
        compare_predictions(nussinov.fold("GGGGAAAACCCC"), thermodynamic.fold("GGGGAAAACCCA"))


# -------------------------------- export -------------------------------- #
def test_fasta_wrapping():
    text = fasta("test", "A" * 130)
    lines = text.strip().split("\n")
    assert lines[0] == ">test" and [len(l) for l in lines[1:]] == [60, 60, 10]


def test_dot_bracket_text():
    seq = "GGGGAAAACCCC"
    assert dot_bracket_text("x", thermodynamic.fold(seq)).splitlines()[2] == "((((....)))) ( -4.18)"
    assert dot_bracket_text("x", nussinov.fold(seq)).splitlines()[2] == "((((....)))) [4 pairs]"


def test_pairs_csv():
    text = pairs_csv(thermodynamic.fold("GGGGAAAACCCC"))
    lines = text.strip().split("\n")
    assert lines[0] == "Position i,Base i,Position j,Base j,Pair,Type"
    assert lines[1] == "1,G,12,C,G-C,Watson-Crick"
    assert len(lines) == 5
    assert pairs_csv(thermodynamic.fold("AAAAAAA")).count("\n") == 1   # header only


def test_mutation_and_scan_csv():
    a = analyse_mutation("GGGGAAAACCCC", 12, "A")
    text = mutation_csv(a)
    assert "Lost pair,1-12,G-C" in text and "C12A" in text
    rows = scan_mutations("GGGAAACCC")
    assert scan_csv(rows).count("\n") == len(rows) + 1


def test_aligned_block_wraps_with_numbers():
    block = aligned_block({"Seq": "A" * 70, "Str": "." * 70}, width=60)
    lines = block.split("\n")
    assert lines[0].split() == ["1", "11", "21", "31", "41", "51"]
    assert lines[1] == "Seq  " + "A" * 60
    assert lines[4].split() == ["61"]
    assert lines[5] == "Seq  " + "A" * 10


def test_summary_report_contains_everything():
    seq = "GGGGAAAACCCC"
    results = {"Nussinov": nussinov.fold(seq), "Thermodynamic": thermodynamic.fold(seq)}
    comp = compare_predictions(results["Nussinov"], results["Thermodynamic"])
    a = analyse_mutation(seq, 12, "A")
    report = summary_report("demo", seq, results, a, comp, scan_mutations(seq))
    for text in ("Length: 12 nt", "-4.18 kcal/mol", "Mutation analysis: C12A", "+3.26",
                 "Mutation scan", "Disclaimer"):
        assert text in report


def test_bundle_zip_and_safe_name():
    data = bundle_zip({"a.txt": "hello", "b.csv": "x,y\n"})
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        assert sorted(z.namelist()) == ["a.txt", "b.csv"]
        assert z.read("a.txt") == b"hello"
    assert safe_name("Yeast tRNA-Phe, 76 nt") == "Yeast_tRNA-Phe_76_nt"
    assert len(safe_name("x" * 100)) == 40 and safe_name("???") == "sequence"
