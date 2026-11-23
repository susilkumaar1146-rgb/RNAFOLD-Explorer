"""
ui/components.py
================
Small reusable pieces of the interface: styling, page headers, the sidebar
status box, tables, plots and the "no sequence yet" guard.
"""

import inspect

import pandas as pd
import streamlit as st

from modules.sequence_validation import gc_content
from ui.state import current_sequence

PAGE_SEQUENCE = "app_pages/sequence_analysis.py"
PAGE_PREDICTION = "app_pages/structure_prediction.py"
PAGE_MUTATION = "app_pages/mutation_analysis.py"
PAGE_VISUALIZATION = "app_pages/visualization.py"
PAGE_RESULTS = "app_pages/results_downloads.py"

_CSS = """
<style>
.block-container {padding-top: 2.2rem; max-width: 1250px;}
div[data-testid="stMetric"] {
    background: #ffffff; border: 1px solid rgba(11,11,11,0.10);
    border-radius: 8px; padding: 10px 14px;
}
div[data-testid="stMetricLabel"] p {color: #52514e;}
h1 {letter-spacing: -0.01em;}
.rfe-subtitle {color: #52514e; font-size: 1.05rem; margin-top: -0.6rem; margin-bottom: 1.2rem;}
.rfe-note {color: #52514e; font-size: 0.9rem;}
</style>
"""


def apply_style() -> None:
    st.markdown(_CSS, unsafe_allow_html=True)


def page_header(title: str, subtitle: str = "") -> None:
    st.title(title)
    if subtitle:
        st.markdown(f'<p class="rfe-subtitle">{subtitle}</p>', unsafe_allow_html=True)


def sidebar_status() -> None:
    """Shows which sequence is loaded, on every page."""
    st.markdown("**Current sequence**")
    seq = current_sequence()
    if seq is None:
        st.caption("None loaded yet - start on *Sequence Analysis*.")
        return
    name = st.session_state.get("sequence_name") or "sequence"
    st.caption(f"{name}\n\n{len(seq)} nt · GC {gc_content(seq):.1f}%")
    results = st.session_state.get("results", {})
    done = [a for a, r in results.items() if r.sequence == seq]
    if done:
        st.caption("Predicted: " + ", ".join(done))


def require_sequence() -> str:
    """Return the current sequence, or explain what to do and stop the page."""
    seq = current_sequence()
    if seq is None:
        st.info("No sequence loaded yet. Enter, upload or choose an example sequence first.")
        st.page_link(PAGE_SEQUENCE, label="Go to Sequence Analysis", icon=":material/arrow_forward:")
        st.stop()
    return seq


def show_plot(fig, key: str) -> None:
    """
    Full-width Plotly chart with our own colours (Streamlit theme off).
    Every chart gets a unique key so two similar charts never clash.
    """
    params = inspect.signature(st.plotly_chart).parameters
    if "width" in params:                      # newer Streamlit: width="stretch" is the default
        st.plotly_chart(fig, theme=None, key=key)
    else:                                      # older Streamlit
        st.plotly_chart(fig, use_container_width=True, theme=None, key=key)


def show_table(rows, height=None) -> None:
    df = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
    if df.empty:
        st.caption("Nothing to show.")
        return
    if height is None:
        st.dataframe(df, hide_index=True)
    else:
        st.dataframe(df, hide_index=True, height=height)


def show_block(text: str) -> None:
    """Monospaced block (sequence / dot-bracket alignments)."""
    st.code(text, language=None)


def note(text: str) -> None:
    st.markdown(f'<p class="rfe-note">{text}</p>', unsafe_allow_html=True)
