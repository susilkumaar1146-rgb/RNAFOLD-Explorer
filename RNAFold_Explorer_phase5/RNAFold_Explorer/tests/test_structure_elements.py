"""Tests for modules/structure_elements.py"""

from modules.dot_bracket import dot_bracket_to_pairs
from modules.structure_elements import (
    LABEL_BULGE, LABEL_EXTERIOR, LABEL_HAIRPIN, LABEL_INTERNAL, LABEL_MULTI,
    LABEL_STEM, analyse_structure,
)


def elements(seq, structure):
    return analyse_structure(seq, dot_bracket_to_pairs(structure))


def test_simple_hairpin():
    e = elements("GGGGAAAACCCC", "((((....))))")
    assert e.counts() == {"Stems": 1, "Hairpin loops": 1, "Bulge loops": 0,
                          "Internal loops": 0, "Multiloops": 0}
    assert e.stems[0].length == 4
    assert e.stem_table()[0]["5' strand"] == "1-4"
    assert e.stem_table()[0]["3' strand"] == "9-12"
    hairpin = [l for l in e.loops if l.label == LABEL_HAIRPIN][0]
    assert hairpin.unpaired_positions == [4, 5, 6, 7]
    assert e.position_labels == [LABEL_STEM] * 4 + [LABEL_HAIRPIN] * 4 + [LABEL_STEM] * 4


def test_bulge_splits_stem():
    e = elements("GGGGAGGGGAAAACCCCCCCC", "((((.((((....))))))))")
    c = e.counts()
    assert c["Stems"] == 2 and c["Bulge loops"] == 1 and c["Hairpin loops"] == 1
    assert e.position_labels[4] == LABEL_BULGE
    assert e.position_element[4] == "B1"


def test_internal_loop():
    e = elements("GGGGAAGGGGAAAACCCCAACCCC", "((((..((((....))))..))))")
    assert e.counts()["Internal loops"] == 1
    internal = [l for l in e.loops if l.label == LABEL_INTERNAL][0]
    assert internal.unpaired_positions == [4, 5, 18, 19]
    assert e.loop_table()[0]["Unpaired positions"] == "5-6, 19-20"


def test_multiloop_and_exterior():
    structure = "..((((..((((....)))).((((((....))))))..))))..."
    seq = "A" * len(structure)        # sequence content not needed for element finding
    e = elements(seq, structure)
    c = e.counts()
    assert c["Stems"] == 3 and c["Multiloops"] == 1 and c["Hairpin loops"] == 2
    exterior = [l for l in e.loops if l.label == LABEL_EXTERIOR][0]
    assert exterior.unpaired_positions == [0, 1, 43, 44, 45]   # 2 + 41 + 3 = 46 nt
    assert e.position_labels[0] == LABEL_EXTERIOR
    assert e.position_labels[6] == LABEL_MULTI


def test_unstructured():
    e = elements("AAAAAA", "......")
    assert e.counts()["Stems"] == 0
    assert set(e.position_labels) == {LABEL_EXTERIOR}


def test_every_position_labelled():
    structure = "((((.((..((...))..))...))))..((...))"
    e = elements("A" * len(structure), structure)
    assert len(e.position_labels) == len(structure)
    for k, ch in enumerate(structure):
        assert (e.position_labels[k] == LABEL_STEM) == (ch != ".")
