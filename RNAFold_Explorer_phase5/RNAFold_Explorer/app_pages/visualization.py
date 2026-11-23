"""Page 5 - Visualization: arc diagram, 2-D schematic, circle plot, elements, DP matrix."""

import streamlit as st

from modules import visualization as V
from modules.dot_bracket import pair_partner_map
from modules.structure_elements import analyse_structure
from ui.components import page_header, require_sequence, show_plot, show_table
from ui.state import available_algorithms, cached_dp_matrix, get_prediction, persistent

DP_MATRIX_LIMIT = 150

page_header("Visualization", "Interactive views of the predicted secondary structure")
seq = require_sequence()
n = len(seq)

algorithms = available_algorithms(seq)
default_algorithm = "Thermodynamic" if "Thermodynamic" in algorithms else "Nussinov"
algorithm = persistent(st.radio, "Structure to display", "viz_algorithm", default_algorithm,
                       options=algorithms, horizontal=True)
with st.spinner("Folding..."):
    result = get_prediction(algorithm)
elements = analyse_structure(seq, result.pairs)

summary = f"{algorithm} prediction · {result.pair_count} base pairs"
if result.energy is not None:
    summary += f" · estimated ΔG {result.energy:.2f} kcal/mol"
st.caption(summary + " · hover over a nucleotide for details; drag to zoom, double-click to reset.")

tab_arc, tab_2d, tab_circle, tab_el, tab_dp = st.tabs(
    ["Arc diagram", "2-D schematic", "Circle plot", "Structural elements", "DP matrix"])

with tab_arc:
    show_plot(V.arc_diagram(seq, result.pairs, f"Arc diagram ({algorithm})"), key="viz_arc")
    st.caption("Nucleotides lie on the axis in sequence order; each arc joins two paired bases. "
               "Nested arcs form stems. Solid arcs: Watson-Crick pairs; dashed arcs: G-U wobble pairs.")
with tab_2d:
    show_plot(V.planar_plot(seq, result.pairs, height=650), key="viz_planar")
    st.caption("Schematic secondary-structure drawing: stems are drawn as ladders and each loop as a "
               "regular polygon. It is a 2-D diagram of base pairing, not an experimentally "
               "determined or 3-D structure. Very large or highly branched structures may overlap.")
with tab_circle:
    show_plot(V.circle_plot(seq, result.pairs, f"Circle plot ({algorithm})", height=600), key="viz_circle")
    st.caption("Bases are placed clockwise from the 5' end; curves join paired bases. In a nested "
               "structure (no pseudoknots) the curves never cross.")
with tab_el:
    counts = elements.counts()
    cols = st.columns(len(counts))
    for col, (label, value) in zip(cols, counts.items()):
        col.metric(label, value)
    st.markdown("**Stems**")
    show_table(elements.stem_table())
    st.markdown("**Loops**")
    show_table(elements.loop_table())
    with st.expander("Per-nucleotide table"):
        partner = pair_partner_map(result.pairs, n)
        show_table([{"Position": k + 1, "Base": seq[k],
                     "Partner": (partner[k] + 1) if partner[k] >= 0 else "-",
                     "Element": elements.position_labels[k], "Element ID": elements.position_element[k]}
                    for k in range(n)], height=400)
    st.caption("Colours and marker shapes in every figure: stem = blue circle, hairpin loop = orange "
               "square, bulge = aqua diamond, internal loop = yellow triangle, multiloop = magenta "
               "pentagon, exterior (unpaired, outside all stems) = grey circle.")
with tab_dp:
    if n > DP_MATRIX_LIMIT:
        st.info(f"The DP matrix view is shown for sequences up to {DP_MATRIX_LIMIT} nt.")
    else:
        show_plot(V.dp_matrix_heatmap(cached_dp_matrix(seq), seq, height=650), key="viz_dp")
        st.caption("Each cell N(i, j) is the largest number of base pairs possible between positions i "
                   "and j (Nussinov). The table is filled from the diagonal (short segments) outwards; "
                   "the top-right cell is the answer for the whole sequence. This view always shows "
                   "the Nussinov table, whichever structure is selected above.")
