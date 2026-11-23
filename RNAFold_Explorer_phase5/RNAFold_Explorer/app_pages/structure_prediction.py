"""Page 3 - Structure Prediction: Nussinov, thermodynamic model and their comparison."""

import streamlit as st

from modules import thermodynamic
from modules import visualization as V
from modules.comparison import compare_predictions
from modules.energy_eval import evaluate_structure
from modules.results_export import aligned_block
from modules.structure_elements import analyse_structure
from ui.components import PAGE_MUTATION, PAGE_VISUALIZATION, page_header, require_sequence, show_block, show_plot, show_table
from ui.state import algorithm_available, get_prediction, persistent

page_header("Structure Prediction", "Predict the secondary structure with one or both algorithms")
seq = require_sequence()

CHOICES = {"Both algorithms": ["Nussinov", "Thermodynamic"],
           "Nussinov only": ["Nussinov"],
           "Thermodynamic model only": ["Thermodynamic"]}
options = list(CHOICES)
if not algorithm_available("Thermodynamic", seq):
    options = ["Nussinov only"]
    st.warning(f"This sequence ({len(seq)} nt) is longer than the thermodynamic model's limit of "
               f"{thermodynamic.MAX_LENGTH} nt, so only Nussinov is available.")

choice = persistent(st.radio, "Algorithm", "predict_algorithms", options[0], options=options,
                    horizontal=True)
if st.button("Run prediction", type="primary", key="btn_predict"):
    st.session_state.predict_choice = CHOICES[choice]

selected = st.session_state.predict_choice
if not selected:
    st.info("Choose the algorithm(s) and press **Run prediction**.")
    st.stop()
selected = [a for a in selected if algorithm_available(a, seq)]

results = {}
with st.spinner("Folding the sequence..."):
    for algorithm in selected:
        results[algorithm] = get_prediction(algorithm)


def render(result) -> None:
    elements = analyse_structure(seq, result.pairs)
    counts = elements.counts()
    c = st.columns(4)
    c[0].metric("Base pairs", result.pair_count)
    if result.energy is not None:
        c[1].metric("Estimated ΔG", f"{result.energy:.2f} kcal/mol",
                    help="Minimum estimated free energy from the simplified model (37 °C).")
    else:
        count = result.optimal_structure_count
        c[1].metric("Equally optimal structures", f"{count:,}" if count else "not counted",
                    help="Number of different structures with the same maximum number of pairs.")
    c[2].metric("Stems", counts["Stems"])
    c[3].metric("Loops (hairpin / bulge / internal / multi)",
                f"{counts['Hairpin loops']} / {counts['Bulge loops']} / "
                f"{counts['Internal loops']} / {counts['Multiloops']}")

    st.markdown("**Dot-bracket structure**")
    show_block(aligned_block({"Sequence": seq, "Structure": result.dot_bracket}))
    for text in result.notes:
        st.caption(text)

    left, right = st.columns(2)
    with left:
        st.markdown(f"**Base pairs ({result.pair_count})**")
        show_table(result.pair_table(), height=320)
    with right:
        st.markdown("**Stems**")
        show_table(elements.stem_table())
        st.markdown("**Loops**")
        show_table(elements.loop_table())

    key = result.algorithm.lower()
    if result.energy is not None:
        st.markdown("**Where the energy comes from**")
        show_plot(V.loop_energy_bar(seq, result.pairs), key=f"energy_bar_{key}")
    else:
        with st.expander("Score this Nussinov structure with the energy model"):
            total = evaluate_structure(seq, result.pairs).total
            st.markdown(f"Under the simplified thermodynamic model this structure has an estimated "
                        f"energy of **{total:.2f} kcal/mol**. This is a separate evaluation for "
                        "comparison; the Nussinov algorithm itself does not use energies.")
            show_plot(V.loop_energy_bar(seq, result.pairs), key=f"energy_bar_{key}")


names = list(results)
tabs = st.tabs(names + (["Comparison"] if len(results) == 2 else []))
for tab, algorithm in zip(tabs, names):
    with tab:
        render(results[algorithm])

if len(results) == 2:
    with tabs[-1]:
        comp = compare_predictions(results["Nussinov"], results["Thermodynamic"])
        nus_partner = {**{i: j for i, j in comp.nussinov.pairs}, **{j: i for i, j in comp.nussinov.pairs}}
        th_partner = {**{i: j for i, j in comp.thermodynamic.pairs},
                      **{j: i for i, j in comp.thermodynamic.pairs}}
        agreement = "".join("|" if nus_partner.get(k) == th_partner.get(k) else "x" for k in range(len(seq)))
        c = st.columns(4)
        c[0].metric("Shared pairs", len(comp.shared))
        c[1].metric("Only in Nussinov", len(comp.only_nussinov))
        c[2].metric("Only in thermodynamic", len(comp.only_thermodynamic))
        c[3].metric("Positions in agreement", f"{comp.position_agreement:.0f} %")
        show_block(aligned_block({"Sequence": seq, "Nussinov": comp.nussinov.dot_bracket,
                                  "Thermo": comp.thermodynamic.dot_bracket, "Agree": agreement}))
        st.caption("Agree: | = same pairing partner (or unpaired in both), x = different.")
        show_table(comp.table_rows())
        st.markdown("**Why the predictions can differ**")
        for line in comp.explanation():
            st.markdown(f"- {line}")
        a1, a2 = st.columns(2)
        with a1:
            show_plot(V.arc_diagram(seq, comp.nussinov.pairs, f"Nussinov ({comp.nussinov.pair_count} pairs)",
                                    height=340), key="cmp_arc_nussinov")
        with a2:
            show_plot(V.arc_diagram(seq, comp.thermodynamic.pairs,
                                    f"Thermodynamic ({comp.thermodynamic.energy:.2f} kcal/mol)", height=340),
                      key="cmp_arc_thermo")

with st.expander("How the two algorithms work"):
    st.markdown(
        "**Nussinov** fills a table N(i, j) with the largest number of pairs possible between "
        "positions i and j, from short segments to long ones, then traces back which choice "
        "(i unpaired, j unpaired, i pairs with j, or a split into two parts) produced the best "
        "value.\n\n"
        "**Thermodynamic model** fills three tables: V(i, j) = lowest energy if i pairs with j "
        "(closing a hairpin, a stack, a bulge/internal loop of at most 30 nt, or a multiloop); "
        "WM(i, j) = lowest energy of a region inside a multiloop; W(j) = lowest energy of the "
        "first j bases. Energies come from a subset of the Turner 2004 nearest-neighbour "
        "parameters (see data/energy_parameters.json). After the traceback, the energy of the "
        "chosen structure is recalculated loop by loop as a correctness check.")

n1, n2 = st.columns(2)
with n1:
    st.page_link(PAGE_VISUALIZATION, label="Explore the structure graphically", icon=":material/insights:")
with n2:
    st.page_link(PAGE_MUTATION, label="Analyse a mutation", icon=":material/compare_arrows:")
