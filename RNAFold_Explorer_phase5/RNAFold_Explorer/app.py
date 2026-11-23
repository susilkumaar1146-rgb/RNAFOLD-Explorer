"""
app.py - RNAFold Explorer (Streamlit entry point)
================================================
Start the web application from the RNAFold_Explorer folder:

    streamlit run app.py

A browser tab opens at http://localhost:8501. Press Ctrl+C in the terminal
to stop the server.

This file only sets up the page, the shared session state and the
navigation menu; each page lives in app_pages/.
"""

import streamlit as st

st.set_page_config(page_title="RNAFold Explorer", page_icon="🧬", layout="wide",
                   initial_sidebar_state="expanded")

from ui.components import apply_style, sidebar_status  # noqa: E402
from ui.state import init_state                        # noqa: E402

init_state()
apply_style()

pages = [
    st.Page("app_pages/home.py", title="Home", icon=":material/home:", default=True),
    st.Page("app_pages/sequence_analysis.py", title="Sequence Analysis", icon=":material/genetics:"),
    st.Page("app_pages/structure_prediction.py", title="Structure Prediction", icon=":material/account_tree:"),
    st.Page("app_pages/mutation_analysis.py", title="Mutation Analysis", icon=":material/compare_arrows:"),
    st.Page("app_pages/visualization.py", title="Visualization", icon=":material/insights:"),
    st.Page("app_pages/results_downloads.py", title="Results & Downloads", icon=":material/download:"),
]
navigation = st.navigation(pages)

with st.sidebar:
    st.divider()
    sidebar_status()

navigation.run()
