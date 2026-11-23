"""Page 2 - Sequence Analysis: input (text / FASTA / examples), validation, composition."""

import os

import streamlit as st

from modules import nussinov, thermodynamic
from modules.mutation_analysis import SCAN_MAX_LENGTH
from modules.results_export import aligned_block
from modules.sequence_validation import base_composition, gc_content, parse_fasta, validate_sequence
from ui.components import PAGE_PREDICTION, page_header, show_block
from ui.state import current_sequence, flash, persistent, set_sequence, show_flash

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SAMPLE_FILE = os.path.join(HERE, "data", "sample_sequences.fasta")


@st.cache_data(show_spinner=False)
def load_examples():
    with open(SAMPLE_FILE, encoding="utf-8") as handle:
        return parse_fasta(handle.read())


def accept(raw: str, name: str, source: str) -> None:
    """Validate; on success store the sequence and re-run the page."""
    check = validate_sequence(raw, max_length=nussinov.MAX_LENGTH)
    if not check.is_valid:
        for error in check.errors:
            st.error(error)
        return
    set_sequence(check.sequence, name.strip() or "my_sequence", source)
    message = f"Sequence **{name}** loaded ({check.length} nt)."
    if check.warnings:
        message += " " + " ".join(check.warnings)
    flash(message)
    st.rerun()


page_header("Sequence Analysis", "Enter an RNA sequence, check it and look at its composition")
show_flash()

method = persistent(st.radio, "Input method", "input_method", "Type or paste",
                    options=["Type or paste", "Upload FASTA file", "Example sequences"],
                    horizontal=True)

with st.container(border=True):
    if method == "Type or paste":
        name = persistent(st.text_input, "Sequence name", "input_name", "my_sequence")
        raw = persistent(st.text_area, "RNA sequence (A, U, G, C - T is converted to U; spaces and "
                         "line breaks are ignored)", "input_raw", "", height=150,
                         placeholder="e.g. GGGGAAAACCCC")
        if st.button("Validate and use this sequence", type="primary", key="btn_use_text"):
            accept(raw, name, "typed")

    elif method == "Upload FASTA file":
        uploaded = st.file_uploader("FASTA file (.fasta, .fa, .fas, .fna or .txt)",
                                    type=["fasta", "fa", "fas", "fna", "txt"], key="input_file")
        if uploaded is not None:
            text = uploaded.getvalue().decode("utf-8", errors="replace")
            records = parse_fasta(text)
            if not records:
                st.error("No sequences were found in this file.")
            else:
                labels = [f"{i + 1}. {header} ({len(''.join(seq.split()))} characters)"
                          for i, (header, seq) in enumerate(records)]
                choice = st.selectbox(f"{len(records)} record(s) found - choose one", labels,
                                      key="input_record")
                header, raw = records[labels.index(choice)]
                if st.button("Validate and use this record", type="primary", key="btn_use_fasta"):
                    accept(raw, header, f"FASTA file {uploaded.name}")

    else:
        examples = load_examples()
        names = [header for header, _ in examples]
        choice = persistent(st.selectbox, "Example sequence", "input_example", names[0], options=names)
        raw = dict(examples)[choice]
        show_block(aligned_block({"Sequence": raw}))
        st.caption("Synthetic examples were designed for this project to show one structural "
                   "feature each; the tRNA is a natural sequence written with unmodified bases.")
        if st.button("Use this example", type="primary", key="btn_use_example"):
            accept(raw, choice, "example")

seq = current_sequence()
st.subheader("Current sequence")
if seq is None:
    st.info("No sequence loaded yet.")
    st.stop()

name = st.session_state.sequence_name
st.markdown(f"**{name}**  ·  source: {st.session_state.sequence_source}")
comp = base_composition(seq)
cols = st.columns(6)
cols[0].metric("Length", f"{len(seq)} nt")
cols[1].metric("GC content", f"{gc_content(seq):.1f} %")
for col, base in zip(cols[2:], "AUGC"):
    col.metric(base, comp[base], help=f"{100 * comp[base] / len(seq):.1f} % of the sequence")
show_block(aligned_block({"5'": seq}))

notes = []
if len(seq) > thermodynamic.MAX_LENGTH:
    notes.append(("warning", f"The thermodynamic model is limited to {thermodynamic.MAX_LENGTH} nt; "
                             "only the Nussinov algorithm can be used for this sequence."))
elif len(seq) > thermodynamic.WARN_LENGTH:
    notes.append(("info", f"Sequences longer than {thermodynamic.WARN_LENGTH} nt can take several "
                          "seconds to fold with the thermodynamic model."))
if len(seq) > SCAN_MAX_LENGTH:
    notes.append(("info", f"The mutation scan (all single mutants) is available up to "
                          f"{SCAN_MAX_LENGTH} nt; single mutations can still be analysed."))
for kind, text in notes:
    (st.warning if kind == "warning" else st.info)(text)

with st.expander("How the input is checked"):
    st.markdown(
        "* Letters are converted to upper case; spaces, tabs and line breaks are removed.\n"
        "* **T is converted to U** so DNA-style input can be used (you are told how many).\n"
        "* Any other character (N, numbers, gaps such as '-') is **rejected with its position**. "
        "Characters are never silently removed, because that would shift every base-pair "
        "position.\n"
        f"* At least 5 nt are required (one pair enclosing a 3-nt hairpin); the maximum is "
        f"{nussinov.MAX_LENGTH} nt.")

st.page_link(PAGE_PREDICTION, label="Next: predict the secondary structure", icon=":material/arrow_forward:")
