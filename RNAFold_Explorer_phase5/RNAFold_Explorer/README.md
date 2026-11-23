# RNAFold Explorer

Interactive RNA secondary structure prediction, visualisation and mutation
impact analysis in Python.

Final-year B.Tech Biotechnology project — Sathvika A S, Susil Kumaar N.

> **Development status:** Phases 2–5 of 7 complete: core algorithms,
> structural elements, mutation analysis, interactive figures and the
> Streamlit web application (123 tests). Phase 6 (testing, benchmarking,
> debugging) and Phase 7 (final documentation) follow.

## Starting the web application

```bat
cd RNAFold_Explorer
.venv\Scripts\activate
streamlit run app.py
```
Your browser opens at http://localhost:8501 (if it does not, open that address
yourself). Press **Ctrl + C** in the Command Prompt to stop the app.

| Page | What it does |
|---|---|
| Home | Project overview, RNA structure background, the two algorithms, limitations |
| Sequence Analysis | Type/paste, upload FASTA or pick an example; validation; length, GC and base counts |
| Structure Prediction | Run Nussinov, the thermodynamic model or both; dot-bracket, pair tables, stems/loops, loop energies, comparison |
| Mutation Analysis | Choose a position and new base; WT vs mutant pairs, ΔΔG, interpretation, figures; all-mutant scan |
| Visualization | Arc diagram, 2-D schematic, circle plot, structural elements, Nussinov DP matrix |
| Results & Downloads | FASTA, dot-bracket, CSV tables, mutation results, text report, everything as ZIP |

## What works now

| Module | Purpose |
|---|---|
| `modules/sequence_validation.py` | Cleans input (upper-case, whitespace, T→U), rejects invalid characters with their positions, parses FASTA, GC content |
| `modules/pairing.py` | Base-pairing rules (A-U, G-C, G-U) and the 3-nt minimum hairpin loop, shared by every algorithm |
| `modules/dot_bracket.py` | Pairs ↔ dot-bracket conversion and structure validation (nesting, duplicates, pseudoknots) |
| `modules/nussinov.py` | Nussinov maximum-base-pair DP, traceback, exact count of co-optimal structures |
| `modules/energy_parameters.py` | Loads `data/energy_parameters.json`; one function per loop type |
| `modules/thermodynamic.py` | Simplified Zuker-style minimum-free-energy DP (V, WM, W states) with traceback |
| `modules/energy_eval.py` | Independent loop-by-loop energy of any structure (checks the DP; scores Nussinov structures) |
| `modules/fold_result.py` | Common result container returned by both algorithms |
| `modules/structure_elements.py` | Finds stems, hairpin/bulge/internal/multi/exterior loops and labels every nucleotide (Phase 3) |
| `modules/structure_layout.py` | Coordinates for the arc diagram, circle plot and the schematic 2-D (loop-polygon) layout (Phase 4) |
| `modules/visualization.py` | Plotly figures: arc diagram, circle plot, 2-D schematic, WT-vs-mutant arcs and drawings, mutation-scan heatmap, loop-energy bars, Nussinov DP matrix (Phase 4) |
| `modules/comparison.py` | Nussinov vs thermodynamic comparison table, shared/unique pairs, explanation (Phase 5) |
| `modules/results_export.py` | FASTA, dot-bracket, CSV, text report and ZIP builders; aligned text blocks (Phase 5) |
| `app.py`, `app_pages/`, `ui/` | Streamlit interface: navigation, six pages, session state, layout helpers (Phase 5) |
| `modules/mutation_analysis.py` | Creates single-nucleotide mutants, compares WT vs mutant (retained/lost/new pairs), ΔΔG, rule-based interpretation, all-mutant scan (Phase 3) |

## Installation (Windows)

```bat
cd RNAFold_Explorer
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

## Running the tests

```bat
python -m pytest -v
```

## Command-line demo

```bat
python demo_cli.py                  :: folds every sequence in data\sample_sequences.fasta
python demo_cli.py GGGGAAAACCCC     :: folds your own sequence
python demo_cli.py --mutate GGGGAAAACCCC C12A   :: mutation analysis
```

## Interactive figures (no Streamlit needed)

```bat
python make_figures.py                       :: example sequence + mutation C10A
python make_figures.py GGGGAAAACCCC C12A     :: your own sequence and mutation
```
The HTML files are written to the `figures` folder; open them in a browser
(an internet connection is needed to load the Plotly library).

ΔΔG = G(mutant) − G(wild type). Positive means the mutant has a higher
(less negative) estimated folding free energy, i.e. it is predicted to be less
stable under this model; this alone says nothing about biological function.

## Energy model (summary)

The model uses a subset of the Turner 2004 nearest-neighbour parameters:
Watson-Crick and G-U stacking, terminal AU/GU penalty, and hairpin, bulge
and internal-loop initiation by size, plus the internal-loop asymmetry and
AU/GU closure terms. Every value in `data/energy_parameters.json` names its
source and is marked `published`, `published_prediction_value` or
`approximate`.

Approximate (non-Turner) values:
* 1×1 and 1×2 internal loops use the n = 4 initiation value (Turner uses lookup tables).
* Multiloop: linear model a = 3.4, b = 0, c = 0.4.

Not modelled: terminal mismatches, dangling ends, special hairpin loops,
special C bulges, small internal-loop tables, coaxial stacking, the
5'GGUC/3'CUGG context, and pseudoknots. All energies are **computational
estimates** and will differ from ViennaRNA/RNAstructure.

## Limits

| Algorithm | Complexity | Warning above | Hard limit |
|---|---|---|---|
| Nussinov | O(n³) | 300 nt | 1000 nt |
| Thermodynamic | O(n³ + n²·30²) | 200 nt | 400 nt |
| Mutation scan (3n folds) | 3n × thermodynamic | — | 80 nt |

## References

1. Nussinov R, Jacobson AB (1980) PNAS 77:6309–6313.
2. Zuker M, Stiegler P (1981) Nucleic Acids Res 9:133–148.
3. Xia T et al. (1998) Biochemistry 37:14719–14735.
4. Mathews DH et al. (2004) PNAS 101:7287–7292.
5. NNDB, Turner 2004: https://rna.urmc.rochester.edu/NNDB/
