"""Page 4 - Mutation Analysis: single-nucleotide mutation and the all-mutant scan."""

import pandas as pd
import streamlit as st

from modules import visualization as V
from modules.dot_bracket import pair_partner_map
from modules.mutation_analysis import interpret
from modules.results_export import aligned_block, scan_csv
from modules.structure_elements import analyse_structure
from ui.components import PAGE_RESULTS, page_header, require_sequence, show_block, show_plot, show_table
from ui.state import (
    SCAN_LIMIT, available_algorithms, cached_mutation, cached_scan, current_mutation,
    current_scan, get_prediction, persistent,
)

page_header("Mutation Analysis",
            "Introduce a single-nucleotide mutation and compare wild-type and mutant predictions")
seq = require_sequence()
n = len(seq)

algorithms = available_algorithms(seq)
default_algorithm = "Thermodynamic" if "Thermodynamic" in algorithms else "Nussinov"
algorithm = persistent(st.radio, "Prediction algorithm (used for both wild type and mutant)",
                       "mut_algorithm", default_algorithm, options=algorithms, horizontal=True)
wild_type = get_prediction(algorithm)
partner = pair_partner_map(wild_type.pairs, n)
elements = analyse_structure(seq, wild_type.pairs)

# ------------------------------ choose mutation ------------------------------ #
with st.container(border=True):
    c1, c2, c3 = st.columns([1, 1, 2])
    with c1:
        position = persistent(st.number_input, "Position (1-based)", "mut_position", 1,
                              min_value=1, max_value=n, step=1)
    original = seq[position - 1]
    with c2:
        replacements = [b for b in "AUGC" if b != original]
        new_base = persistent(st.selectbox, f"Replace {original}{position} with", "mut_base",
                              replacements[0], options=replacements)
    with c3:
        if partner[position - 1] >= 0:
            p = partner[position - 1]
            status = f"paired with {seq[p]}{p + 1}"
        else:
            status = "unpaired"
        st.markdown(f"**Original nucleotide:** {original}{position}  \n"
                    f"In the wild-type {algorithm} prediction: **{status}** "
                    f"({elements.position_labels[position - 1].lower()}).")
        st.markdown(f"**Mutation:** `{original}{position}{new_base}`")

    mutant_seq = seq[:position - 1] + new_base + seq[position:]
    marker = "".join("^" if k == position - 1 else " " for k in range(n))
    show_block(aligned_block({"Wild type": seq, "Mutant": mutant_seq, "": marker}))
    run = st.button("Analyse mutation", type="primary", key="btn_mutation")

if run:
    with st.spinner("Folding wild type and mutant..."):
        st.session_state.mutation = cached_mutation(seq, int(position), new_base, algorithm)

analysis = current_mutation()
if analysis is None:
    st.info("Choose a position and a replacement base, then press **Analyse mutation**.")
else:
    m, comp = analysis.mutation, analysis.comparison
    st.subheader(f"Results for {m.label} ({analysis.algorithm})")
    if (m.position, m.replacement, analysis.algorithm) != (position, new_base, algorithm):
        st.caption("Showing the last analysed mutation; press **Analyse mutation** to update.")

    c = st.columns(5)
    c[0].metric("Base pairs (wild type)", analysis.wild_type.pair_count)
    c[1].metric("Base pairs (mutant)", analysis.mutant.pair_count,
                delta=analysis.mutant.pair_count - analysis.wild_type.pair_count)
    c[2].metric("Retained", len(comp.retained))
    c[3].metric("Lost", len(comp.lost))
    c[4].metric("Newly formed", len(comp.new))
    e = st.columns(3)
    e[0].metric("Estimated ΔG wild type", f"{analysis.wt_energy:.2f} kcal/mol")
    e[1].metric("Estimated ΔG mutant", f"{analysis.mut_energy:.2f} kcal/mol")
    e[2].metric("ΔΔG = G(mutant) − G(wild type)", f"{analysis.ddg:+.2f} kcal/mol",
                delta=f"{analysis.ddg:+.2f}", delta_color="inverse",
                help="Positive: higher estimated free energy (less stable under this model). "
                     "Negative: lower (more stable).")
    st.caption(analysis.energy_note)

    changes = "".join("^" if k == m.index else ("*" if k in comp.changed_positions else " ")
                      for k in range(n))
    show_block(aligned_block({"WT seq": analysis.wild_type.sequence, "WT struct": analysis.wild_type.dot_bracket,
                              "Mut seq": analysis.mutant.sequence, "Mut struct": analysis.mutant.dot_bracket,
                              "Changes": changes}))
    st.caption("^ = mutated position, * = pairing partner changed.")

    st.markdown("**Interpretation**")
    for line in interpret(analysis):
        st.markdown(f"- {line}")

    tab_arc, tab_2d, tab_pairs, tab_el = st.tabs(["Arc comparison", "2-D side by side",
                                                  "Base-pair changes", "Structural elements"])
    with tab_arc:
        show_plot(V.mutation_arc_comparison(analysis), key="mut_arc")
        st.caption("Wild-type pairs above the axis, mutant pairs below. Grey = retained, "
                   "red dashed = lost in the mutant, green = newly formed; the mutated base is circled.")
    with tab_2d:
        fig_wt, fig_mut = V.mutation_planar_plots(analysis, height=480)
        left, right = st.columns(2)
        with left:
            show_plot(fig_wt, key="mut_planar_wt")
        with right:
            show_plot(fig_mut, key="mut_planar_mut")
        st.caption("Schematic 2-D drawings (not 3-D models).")
    with tab_pairs:
        def rows(pairs, sequence):
            return [{"Pair": f"{i + 1}-{j + 1}", "Bases": f"{sequence[i]}-{sequence[j]}"} for i, j in pairs]
        p1, p2, p3 = st.columns(3)
        with p1:
            st.markdown(f"**Lost ({len(comp.lost)})**")
            show_table(rows(comp.lost, analysis.wild_type.sequence))
        with p2:
            st.markdown(f"**Newly formed ({len(comp.new)})**")
            show_table(rows(comp.new, analysis.mutant.sequence))
        with p3:
            st.markdown(f"**Retained ({len(comp.retained)})**")
            show_table(rows(comp.retained, analysis.wild_type.sequence), height=300)
    with tab_el:
        show_table([{"Element": k, "Wild type": analysis.wt_elements[k], "Mutant": analysis.mut_elements[k],
                     "Change": analysis.mut_elements[k] - analysis.wt_elements[k]} for k in analysis.wt_elements])
        st.markdown(f"Mutated site in the wild type: {analysis.wt_site}.  \n"
                    f"Mutated site in the mutant: {analysis.mut_site}.")

# ------------------------------ mutation scan ------------------------------ #
st.divider()
st.subheader("Mutation scan")
st.markdown("Folds **every** possible single-nucleotide mutant (3 per position) with the selected "
            "algorithm and maps ΔΔG along the sequence, to show which positions the predicted "
            "structure is most sensitive to.")
if n > SCAN_LIMIT:
    st.info(f"The scan is available for sequences up to {SCAN_LIMIT} nt "
            f"(this one has {n} nt, which would need {3 * n} folds).")
else:
    if st.button(f"Run mutation scan ({3 * n} folds)", key="btn_scan"):
        bar = st.progress(0.0, text="Folding mutants...")
        rows = cached_scan(seq, algorithm, _progress=lambda f: bar.progress(min(f, 1.0), text="Folding mutants..."))
        bar.empty()
        st.session_state.scan = {"sequence": seq, "algorithm": algorithm, "rows": rows}

    scan = current_scan()
    if scan is not None:
        rows = scan["rows"]
        st.caption(f"Scan computed with the {scan['algorithm']} algorithm.")
        show_plot(V.mutation_scan_heatmap(rows, seq), key="scan_heatmap")
        df = pd.DataFrame(rows)
        t1, t2 = st.columns(2)
        with t1:
            st.markdown("**Largest ΔΔG (predicted least stable mutants)**")
            show_table(df.sort_values("ΔΔG (kcal/mol)", ascending=False).head(10)[
                ["Mutation", "ΔΔG (kcal/mol)", "Pairs lost", "Pairs gained"]])
        with t2:
            st.markdown("**Smallest ΔΔG (predicted most stable mutants)**")
            show_table(df.sort_values("ΔΔG (kcal/mol)").head(10)[
                ["Mutation", "ΔΔG (kcal/mol)", "Pairs lost", "Pairs gained"]])
        st.download_button("Download scan results (CSV)", scan_csv(rows),
                           file_name="mutation_scan.csv", mime="text/csv", key="dl_scan_page")

st.page_link(PAGE_RESULTS, label="Download all results", icon=":material/download:")
