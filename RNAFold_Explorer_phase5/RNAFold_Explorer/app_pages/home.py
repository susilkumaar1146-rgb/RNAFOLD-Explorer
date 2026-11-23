"""Page 1 - Home: project overview, biology background and how to use the app."""

import streamlit as st

from ui.components import (
    PAGE_MUTATION, PAGE_PREDICTION, PAGE_RESULTS, PAGE_SEQUENCE, PAGE_VISUALIZATION,
    page_header, show_table,
)

page_header("🧬 RNAFold Explorer",
            "Interactive RNA secondary structure prediction, visualisation and "
            "mutation impact analysis")

st.markdown(
    "RNAFold Explorer predicts how a single RNA strand folds back on itself, using only its "
    "nucleotide sequence. It implements two dynamic programming algorithms in Python, draws "
    "the predicted structures interactively and shows how a single-nucleotide mutation "
    "changes the predicted fold and its estimated stability.")

col1, col2, col3 = st.columns(3)
with col1:
    with st.container(border=True):
        st.markdown("**1 · Predict**")
        st.caption("Nussinov (maximum base pairs) and a simplified thermodynamic model "
                   "(minimum estimated free energy).")
with col2:
    with st.container(border=True):
        st.markdown("**2 · Visualise**")
        st.caption("Arc diagram, circle plot and a schematic 2-D drawing with stems and "
                   "loops coloured.")
with col3:
    with st.container(border=True):
        st.markdown("**3 · Mutate**")
        st.caption("Change one nucleotide and compare wild-type and mutant structures, "
                   "base pairs and ΔΔG.")

st.subheader("Getting started")
st.page_link(PAGE_SEQUENCE, label="1. Enter, upload or choose an RNA sequence", icon=":material/genetics:")
st.page_link(PAGE_PREDICTION, label="2. Predict its secondary structure", icon=":material/account_tree:")
st.page_link(PAGE_VISUALIZATION, label="3. Explore the structure graphically", icon=":material/insights:")
st.page_link(PAGE_MUTATION, label="4. Analyse a single-nucleotide mutation", icon=":material/compare_arrows:")
st.page_link(PAGE_RESULTS, label="5. Download your results", icon=":material/download:")
st.caption("Your sequence and results are kept while you move between pages "
           "(they are cleared when the browser tab is reloaded).")

st.divider()
st.header("RNA secondary structure")
st.markdown(
    "RNA is built from four nucleotides: **adenine (A)**, **uracil (U)**, **guanine (G)** and "
    "**cytosine (C)**. A single RNA strand folds by forming hydrogen-bonded base pairs between "
    "complementary regions of the same molecule. The pattern of these pairs is the "
    "**secondary structure**. The structures of non-coding RNAs such as tRNA, rRNA, miRNA "
    "precursors and riboswitches are closely linked to their biological functions.")

left, right = st.columns(2)
with left:
    st.markdown("**Base pairs considered**")
    show_table([
        {"Pair": "A–U", "Type": "Watson-Crick (canonical)", "Hydrogen bonds": 2},
        {"Pair": "G–C", "Type": "Watson-Crick (canonical)", "Hydrogen bonds": 3},
        {"Pair": "G–U", "Type": "Wobble", "Hydrogen bonds": 2},
    ])
    st.caption("A hairpin loop must contain at least 3 unpaired nucleotides.")
with right:
    st.markdown("**Structural elements**")
    show_table([
        {"Element": "Stem (helix)", "Description": "Consecutive stacked base pairs"},
        {"Element": "Hairpin loop", "Description": "Unpaired bases closed by one pair"},
        {"Element": "Bulge loop", "Description": "Unpaired bases on one side of a stem"},
        {"Element": "Internal loop", "Description": "Unpaired bases on both sides between two stems"},
        {"Element": "Multiloop", "Description": "Junction where three or more stems meet"},
    ])

st.markdown("**Dot-bracket notation** writes a structure as one character per nucleotide: "
            "`(` and `)` mark the two partners of a pair, `.` an unpaired base.")
st.code("Sequence    GGGGAAAACCCC\nStructure   ((((....))))   4 pairs enclosing a 4-nt hairpin loop",
        language=None)
st.caption("Illustrative example. The application always computes structures with the algorithms; "
           "nothing is hard-coded.")

st.divider()
st.header("Prediction algorithms")
a1, a2 = st.columns(2)
with a1:
    with st.container(border=True):
        st.subheader("Nussinov algorithm")
        st.markdown("Finds the nested structure with the **maximum number of base pairs** "
                    "(Nussinov & Jacobson, 1980). Every pair scores +1; energies are not used.")
        st.latex(r"N(i,j)=\max\begin{cases}N(i+1,j)\\N(i,j-1)\\N(i+1,j-1)+1 & \text{if } i\cdot j \text{ can pair}"
                 r"\\\max_{i<k<j} N(i,k)+N(k+1,j)\end{cases}")
        st.caption("Time O(n³), memory O(n²). Traceback recovers the pairs.")
with a2:
    with st.container(border=True):
        st.subheader("Simplified thermodynamic model")
        st.markdown("Finds the structure with the **minimum estimated free energy**, in the style "
                    "of Zuker & Stiegler (1981). Stacked pairs stabilise (negative ΔG); hairpin, "
                    "bulge, internal and multiloops cost energy (positive ΔG).")
        st.latex(r"V(i,j)=\min\{\,H(i,j),\; \Delta G_{\text{loop}}(i,j,p,q)+V(p,q),\; "
                 r"\text{multiloop}\,\}")
        st.caption("Uses a subset of the Turner 2004 nearest-neighbour parameters (37 °C). "
                   "States V, WM and W; time O(n³), memory O(n²).")

st.divider()
st.header("Assumptions and limitations")
st.info(
    "• All predictions are **computational estimates**, not experimentally determined structures.\n\n"
    "• Only **nested** structures are predicted; pseudoknots are not supported.\n\n"
    "• The thermodynamic model is **simplified**: terminal mismatches, dangling ends, special "
    "hairpin loops, small internal-loop tables and coaxial stacking are not modelled, and the "
    "multiloop parameters are approximate. Energies will differ from ViennaRNA or RNAstructure.\n\n"
    "• A change in predicted structure or ΔΔG does not by itself show a loss or gain of "
    "biological function.")

with st.expander("References"):
    st.markdown(
        "1. Nussinov R, Jacobson AB (1980) Fast algorithm for predicting the secondary structure "
        "of single-stranded RNA. *PNAS* 77:6309–6313.\n"
        "2. Zuker M, Stiegler P (1981) Optimal computer folding of large RNA sequences using "
        "thermodynamics and auxiliary information. *Nucleic Acids Res* 9:133–148.\n"
        "3. Xia T et al. (1998) Thermodynamic parameters for an expanded nearest-neighbor model "
        "for formation of RNA duplexes with Watson-Crick base pairs. *Biochemistry* 37:14719–14735.\n"
        "4. Mathews DH et al. (2004) Incorporating chemical modification constraints into a dynamic "
        "programming algorithm for prediction of RNA secondary structure. *PNAS* 101:7287–7292.\n"
        "5. Nearest Neighbor Database (Turner 2004): https://rna.urmc.rochester.edu/NNDB/\n"
        "6. ViennaRNA documentation, dot-bracket notation: "
        "https://viennarna.readthedocs.io/en/latest/io/rna_structures.html")

st.caption("RNAFold Explorer · B.Tech Biotechnology final-year project · "
           "Sathvika A S, Susil Kumaar N")
