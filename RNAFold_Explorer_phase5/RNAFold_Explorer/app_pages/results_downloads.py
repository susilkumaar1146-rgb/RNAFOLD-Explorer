"""Page 6 - Results & Downloads: everything computed so far, as files."""

import streamlit as st

from modules.comparison import compare_predictions
from modules.results_export import (
    bundle_zip, dot_bracket_text, fasta, mutation_csv, pairs_csv, safe_name, scan_csv,
    summary_report,
)
from ui.components import PAGE_MUTATION, PAGE_PREDICTION, page_header, require_sequence, show_table
from ui.state import available_algorithms, current_mutation, current_scan, get_prediction, stored_results

page_header("Results & Downloads", "Download the sequence, predictions and analyses as files")
seq = require_sequence()
name = st.session_state.sequence_name or "sequence"
base = safe_name(name)

results = stored_results()
mutation = current_mutation()
scan = current_scan()

if not results:
    st.info("No structure has been predicted for this sequence yet.")
    if st.button("Run all available predictions now", type="primary", key="btn_results_predict"):
        with st.spinner("Folding..."):
            for algorithm in available_algorithms(seq):
                get_prediction(algorithm)
        st.rerun()
    st.page_link(PAGE_PREDICTION, label="Or go to Structure Prediction", icon=":material/account_tree:")

comparison = None
if "Nussinov" in results and "Thermodynamic" in results:
    comparison = compare_predictions(results["Nussinov"], results["Thermodynamic"])

# ------------------------------ overview ------------------------------ #
st.subheader("Summary")
overview = [{"Item": "Sequence", "Status": f"{name} ({len(seq)} nt)"}]
for algorithm in ("Nussinov", "Thermodynamic"):
    r = results.get(algorithm)
    if r is None:
        overview.append({"Item": f"{algorithm} prediction", "Status": "not run"})
    else:
        text = f"{r.pair_count} pairs"
        if r.energy is not None:
            text += f", {r.energy:.2f} kcal/mol"
        overview.append({"Item": f"{algorithm} prediction", "Status": text})
overview.append({"Item": "Mutation analysis",
                 "Status": f"{mutation.mutation.label} ({mutation.algorithm}), ΔΔG {mutation.ddg:+.2f} kcal/mol"
                 if mutation else "not run"})
overview.append({"Item": "Mutation scan",
                 "Status": f"{len(scan['rows'])} mutants ({scan['algorithm']})" if scan else "not run"})
show_table(overview)
if mutation is None:
    st.page_link(PAGE_MUTATION, label="Run a mutation analysis", icon=":material/compare_arrows:")

# ------------------------------ files ------------------------------ #
files = {f"{base}.fasta": fasta(name, seq)}
for algorithm, r in results.items():
    files[f"{base}_{algorithm.lower()}_structure.txt"] = dot_bracket_text(name, r)
    files[f"{base}_{algorithm.lower()}_pairs.csv"] = pairs_csv(r)
if mutation is not None:
    label = mutation.mutation.label
    files[f"{base}_{label}_mutant.fasta"] = fasta(f"{name} mutant {label}", mutation.mutant.sequence)
    files[f"{base}_{label}_comparison.csv"] = mutation_csv(mutation)
if scan is not None:
    files[f"{base}_mutation_scan.csv"] = scan_csv(scan["rows"])
report = summary_report(name, seq, results, mutation, comparison, scan["rows"] if scan else None)
files[f"{base}_summary_report.txt"] = report

st.subheader("Download files")
st.download_button("Download everything (ZIP)", bundle_zip(files), file_name=f"{base}_RNAFold_Explorer.zip",
                   mime="application/zip", type="primary", key="dl_zip")

DESCRIPTIONS = {
    ".fasta": "FASTA sequence",
    "_structure.txt": "Dot-bracket structure (Vienna-style text)",
    "_pairs.csv": "Base-pair list (CSV)",
    "_comparison.csv": "Wild type vs mutant comparison (CSV)",
    "_mutation_scan.csv": "Mutation scan, all single mutants (CSV)",
    "_summary_report.txt": "Full prediction summary (text)",
}
for k, (filename, content) in enumerate(files.items()):
    description = next((d for suffix, d in DESCRIPTIONS.items() if filename.endswith(suffix)), "File")
    if filename.endswith("_mutant.fasta"):
        description = "Mutant sequence (FASTA)"
    mime = "text/csv" if filename.endswith(".csv") else "text/plain"
    c1, c2 = st.columns([3, 2])
    c1.markdown(f"**{description}**  \n`{filename}`")
    c2.download_button("Download", content, file_name=filename, mime=mime, key=f"dl_{k}_{filename}")

with st.expander("Preview the summary report"):
    st.code(report, language=None)
