"""
Tests for modules/visualization.py.
These need Plotly installed (pip install -r requirements.txt); they are
skipped automatically if Plotly is missing.
"""

import pytest

pytest.importorskip("plotly")

from modules import nussinov, thermodynamic                   # noqa: E402
from modules import visualization as V                         # noqa: E402
from modules.mutation_analysis import analyse_mutation, scan_mutations  # noqa: E402

SEQ = "GGGCCCAAGCUCGAAAGAGCAAUGGACGAAAGUCCAUAAGGGCCC"     # multiloop example
WOBBLE = "GGGUGCAAAAGUACCC"


def names(fig):
    return [t.name for t in fig.data]


def all_hover(fig):
    texts = []
    for t in fig.data:
        if t.hovertext is not None:
            texts.extend(t.hovertext)
    return texts


@pytest.fixture(scope="module")
def fold():
    return thermodynamic.fold(SEQ)


def test_arc_diagram_traces_and_hover(fold):
    fig = V.arc_diagram(SEQ, fold.pairs)
    assert "Watson-Crick pair" in names(fig)
    assert {"Stem", "Hairpin loop", "Multiloop"} <= set(names(fig))
    hover = all_hover(fig)
    assert len(hover) == len(SEQ)                        # one hover per nucleotide
    assert any("Paired with" in h for h in hover) and any("Unpaired" in h for h in hover)


def test_every_nucleotide_drawn_once(fold):
    fig = V.arc_diagram(SEQ, fold.pairs)
    xs = []
    for t in fig.data:
        if t.mode and "markers" in t.mode:
            xs.extend(t.x)
    assert sorted(xs) == list(range(1, len(SEQ) + 1))


def test_wobble_pairs_dashed():
    fig = V.arc_diagram(WOBBLE, thermodynamic.fold(WOBBLE).pairs)
    wobble = [t for t in fig.data if t.name == "G-U wobble pair"][0]
    assert wobble.line.dash == "dash"


def test_circle_and_planar(fold):
    for fig in (V.circle_plot(SEQ, fold.pairs), V.planar_plot(SEQ, fold.pairs)):
        assert len(all_hover(fig)) == len(SEQ)
        fig.to_json()                                    # figure must be serialisable


def test_planar_title_says_schematic(fold):
    fig = V.planar_plot(SEQ, fold.pairs)
    assert "Schematic" in fig.layout.title.text and "not a 3-D" in fig.layout.title.text


def test_mutation_figures():
    a = analyse_mutation(SEQ, 10, "A")
    fig = V.mutation_arc_comparison(a)
    assert "Lost in mutant" in names(fig)
    assert f"Mutation {a.mutation.label}" in names(fig)
    fig_wt, fig_mut = V.mutation_planar_plots(a)
    assert "Mutation site" in names(fig_wt) and "Mutated nucleotide" in names(fig_mut)
    fig.to_json(); fig_wt.to_json(); fig_mut.to_json()


def test_scan_heatmap_and_energy_bar(fold):
    rows = scan_mutations(SEQ[:20])
    fig = V.mutation_scan_heatmap(rows, SEQ[:20])
    assert len(fig.data[0].z) == 4 and len(fig.data[0].z[0]) == 20
    bar = V.loop_energy_bar(SEQ, fold.pairs)
    assert {"Stabilising (ΔG < 0)", "Destabilising (ΔG ≥ 0)"} <= set(names(bar))
    bar.to_json(); fig.to_json()


def test_dp_matrix_heatmap():
    seq = "GGGGAAAACCCC"
    fig = V.dp_matrix_heatmap(nussinov.fill_matrix(seq), seq)
    z = fig.data[0].z
    assert z[0][11] == 4 and z[5][2] is None             # lower triangle left empty
