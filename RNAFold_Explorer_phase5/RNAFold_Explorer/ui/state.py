"""
ui/state.py
===========
Session state and cached computations.

Streamlit re-runs the page script from top to bottom after every click.
Two tools keep the app fast and stop results from disappearing:

* st.session_state  - a dictionary that survives re-runs and page changes.
  We keep the current sequence and the latest results in it.
* st.cache_data      - remembers the output of a function for given inputs,
  so folding the same sequence twice costs nothing.

Keys used in st.session_state
-----------------------------
sequence, sequence_name, sequence_source   the validated sequence in use
results          {"Nussinov": FoldResult, "Thermodynamic": FoldResult}
predict_choice   algorithms the user asked for on the prediction page
mutation         latest MutationAnalysis
scan             {"algorithm": ..., "rows": [...]} from the mutation scan
flash            a message shown once after a page re-run
_saved_<key>     remembered widget values (see persistent())
"""

from typing import List, Optional

import streamlit as st

from modules import nussinov, thermodynamic
from modules.fold_result import FoldResult
from modules.mutation_analysis import SCAN_MAX_LENGTH, analyse_mutation, scan_mutations

DEFAULTS = {
    "sequence": None,
    "sequence_name": None,
    "sequence_source": None,
    "results": {},
    "predict_choice": None,
    "mutation": None,
    "scan": None,
    "flash": None,
}

ALGORITHMS = ["Nussinov", "Thermodynamic"]


def init_state() -> None:
    for key, value in DEFAULTS.items():
        if key not in st.session_state:
            st.session_state[key] = dict(value) if isinstance(value, dict) else value


# ------------------------------------------------------------------ #
# Current sequence
# ------------------------------------------------------------------ #
def set_sequence(sequence: str, name: str, source: str) -> None:
    """Store a validated sequence; results for a different sequence are cleared."""
    if sequence != st.session_state.sequence:
        st.session_state.results = {}
        st.session_state.predict_choice = None
        st.session_state.mutation = None
        st.session_state.scan = None
    st.session_state.sequence = sequence
    st.session_state.sequence_name = name
    st.session_state.sequence_source = source


def current_sequence() -> Optional[str]:
    return st.session_state.get("sequence")


def algorithm_available(algorithm: str, sequence: str) -> bool:
    limit = nussinov.MAX_LENGTH if algorithm == "Nussinov" else thermodynamic.MAX_LENGTH
    return len(sequence) <= limit


def available_algorithms(sequence: str) -> List[str]:
    return [a for a in ALGORITHMS if algorithm_available(a, sequence)]


# ------------------------------------------------------------------ #
# Cached computations
# ------------------------------------------------------------------ #
@st.cache_data(show_spinner=False, max_entries=100)
def cached_fold(sequence: str, algorithm: str) -> FoldResult:
    if algorithm == "Nussinov":
        return nussinov.fold(sequence)
    return thermodynamic.fold(sequence)


@st.cache_data(show_spinner=False, max_entries=100)
def cached_mutation(sequence: str, position: int, new_base: str, algorithm: str):
    wild_type = cached_fold(sequence, algorithm)
    return analyse_mutation(sequence, position, new_base, algorithm, wild_type=wild_type)


@st.cache_data(show_spinner=False, max_entries=20)
def cached_scan(sequence: str, algorithm: str, _progress=None) -> list:
    # The leading underscore tells Streamlit not to hash the progress callback.
    return scan_mutations(sequence, algorithm, progress=_progress)


@st.cache_data(show_spinner=False, max_entries=10)
def cached_dp_matrix(sequence: str):
    return nussinov.fill_matrix(sequence)


def get_prediction(algorithm: str) -> FoldResult:
    """Prediction for the current sequence (computed once, then reused)."""
    sequence = current_sequence()
    result = st.session_state.results.get(algorithm)
    if result is None or result.sequence != sequence:
        result = cached_fold(sequence, algorithm)
        st.session_state.results[algorithm] = result
    return result


def stored_results() -> dict:
    """Predictions already computed for the current sequence."""
    sequence = current_sequence()
    return {a: r for a, r in st.session_state.results.items() if r.sequence == sequence}


def current_mutation():
    a = st.session_state.get("mutation")
    if a is not None and a.wild_type.sequence == current_sequence():
        return a
    return None


def current_scan():
    s = st.session_state.get("scan")
    if s is not None and s.get("sequence") == current_sequence():
        return s
    return None


SCAN_LIMIT = SCAN_MAX_LENGTH


# ------------------------------------------------------------------ #
# Remembered widget values
# ------------------------------------------------------------------ #
def persistent(widget, label: str, key: str, default, options: Optional[list] = None, **kwargs):
    """
    Create a widget whose value is remembered when the user switches pages.
    (Streamlit forgets the state of widgets that are not on the current page.)
    `options` is passed to radio/selectbox; a remembered value that is no longer
    a valid option falls back to `default`.
    """
    store = "_saved_" + key
    if key not in st.session_state:
        value = st.session_state.get(store, default)
        if options is not None and value not in options:
            value = default
        if "min_value" in kwargs and value is not None:
            value = max(kwargs["min_value"], min(kwargs.get("max_value", value), value))
        st.session_state[key] = value
    elif options is not None and st.session_state[key] not in options:
        st.session_state[key] = default
    if options is not None:
        value = widget(label, options, key=key, **kwargs)
    else:
        if "min_value" in kwargs:
            current = st.session_state[key]
            clamped = max(kwargs["min_value"], min(kwargs.get("max_value", current), current))
            if clamped != current:
                st.session_state[key] = clamped
        value = widget(label, key=key, **kwargs)
    st.session_state[store] = value
    return value


# ------------------------------------------------------------------ #
# One-time messages (shown after st.rerun())
# ------------------------------------------------------------------ #
def flash(message: str, kind: str = "success") -> None:
    st.session_state.flash = (kind, message)


def show_flash() -> None:
    item = st.session_state.get("flash")
    if item:
        kind, message = item
        {"success": st.success, "info": st.info, "warning": st.warning}.get(kind, st.info)(message)
        st.session_state.flash = None
