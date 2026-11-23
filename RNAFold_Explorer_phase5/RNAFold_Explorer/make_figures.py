"""
make_figures.py
===============
Saves every Phase 4 figure as an interactive HTML page (open them in any
web browser - zoom, pan and hover work). No Streamlit needed.

Usage (from the RNAFold_Explorer folder):
    python make_figures.py                         # uses the multiloop example
    python make_figures.py GGGGAAAACCCC C12A       # your sequence + a mutation
Output goes to the "figures" folder.
"""

import os
import sys

from modules import nussinov, thermodynamic
from modules import visualization as V
from modules.mutation_analysis import SCAN_MAX_LENGTH, analyse_mutation, parse_mutation, scan_mutations
from modules.sequence_validation import validate_sequence

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
DEFAULT_SEQ = "GGGCCCAAGCUCGAAAGAGCAAUGGACGAAAGUCCAUAAGGGCCC"
DEFAULT_MUT = "C10A"


def main() -> None:
    raw = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SEQ
    mut_text = sys.argv[2] if len(sys.argv) > 2 else (DEFAULT_MUT if len(sys.argv) == 1 else None)
    check = validate_sequence(raw, max_length=thermodynamic.MAX_LENGTH)
    if not check.is_valid:
        print("\n".join("ERROR: " + e for e in check.errors))
        return
    seq = check.sequence
    os.makedirs(OUT, exist_ok=True)

    thermo = thermodynamic.fold(seq)
    nus = nussinov.fold(seq)
    figures = {
        "arc_thermodynamic": V.arc_diagram(seq, thermo.pairs, f"Arc diagram - thermodynamic ({thermo.energy:.2f} kcal/mol)"),
        "arc_nussinov": V.arc_diagram(seq, nus.pairs, f"Arc diagram - Nussinov ({nus.pair_count} pairs)"),
        "circle_thermodynamic": V.circle_plot(seq, thermo.pairs),
        "planar_thermodynamic": V.planar_plot(seq, thermo.pairs),
        "loop_energies": V.loop_energy_bar(seq, thermo.pairs),
    }
    if len(seq) <= 150:
        figures["nussinov_dp_matrix"] = V.dp_matrix_heatmap(nussinov.fill_matrix(seq), seq)
    if mut_text:
        m = parse_mutation(mut_text, seq)
        a = analyse_mutation(seq, m.position, m.replacement, "Thermodynamic", wild_type=thermo)
        figures[f"mutation_{m.label}_arcs"] = V.mutation_arc_comparison(a)
        fig_wt, fig_mut = V.mutation_planar_plots(a)
        figures[f"mutation_{m.label}_planar_wild_type"] = fig_wt
        figures[f"mutation_{m.label}_planar_mutant"] = fig_mut
    if len(seq) <= SCAN_MAX_LENGTH:
        print(f"Running mutation scan ({3 * len(seq)} folds)...")
        figures["mutation_scan"] = V.mutation_scan_heatmap(scan_mutations(seq), seq)

    for name, fig in figures.items():
        path = os.path.join(OUT, name + ".html")
        fig.write_html(path, include_plotlyjs="cdn")
        print("saved", path)


if __name__ == "__main__":
    main()
