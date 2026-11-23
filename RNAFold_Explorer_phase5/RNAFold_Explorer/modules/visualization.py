"""
visualization.py
================
Interactive Plotly figures for RNAFold Explorer. Every function returns a
plotly.graph_objects.Figure, which Streamlit shows with st.plotly_chart().
Zooming, panning and hover tooltips come built into Plotly figures.

Figures
-------
arc_diagram()              arcs over a position axis, nucleotides coloured by element
circle_plot()              bases on a circle, pairs as curved chords
planar_plot()              schematic 2-D secondary-structure drawing
mutation_arc_comparison()  wild-type arcs above the axis, mutant arcs below
mutation_planar_plots()    wild-type and mutant 2-D drawings for side-by-side display
mutation_scan_heatmap()    ddG of every single mutant (position x new base)
loop_energy_bar()          contribution of every loop to the estimated energy
dp_matrix_heatmap()        the Nussinov N(i, j) table (teaching aid)

The coordinates come from structure_layout.py and the element labels from
structure_elements.py, so this module only handles drawing.

Colour choices
--------------
* Structural elements use a fixed categorical order (validated for colour-
  vision deficiency on adjacent slots) and each element also has its own
  marker SHAPE, so identity never relies on colour alone.
* Base-pair lines are neutral ink; G-U wobble pairs are dashed.
* Mutation comparison: retained = grey, lost = red dashed, new = green,
  each named in the legend.
* Energy figures use a diverging blue (stabilising) / red (destabilising)
  scale with a neutral grey midpoint at 0.
"""

import math
from typing import Dict, List, Optional, Sequence, Tuple

import plotly.graph_objects as go

from modules.dot_bracket import pair_partner_map
from modules.energy_eval import evaluate_structure
from modules.pairing import pair_type
from modules.structure_elements import (
    ELEMENT_ORDER, LABEL_BULGE, LABEL_EXTERIOR, LABEL_HAIRPIN, LABEL_INTERNAL,
    LABEL_MULTI, LABEL_STEM, analyse_structure,
)
from modules.structure_layout import arc_points, chord_points, circle_positions, planar_layout

Pair = Tuple[int, int]

# ---------------------------- palette ---------------------------------- #
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
FONT = "system-ui, -apple-system, 'Segoe UI', sans-serif"

ELEMENT_COLORS = {
    LABEL_STEM: "#2a78d6",       # blue
    LABEL_HAIRPIN: "#eb6834",    # orange
    LABEL_BULGE: "#1baf7a",      # aqua
    LABEL_INTERNAL: "#eda100",   # yellow
    LABEL_MULTI: "#e87ba4",      # magenta
    LABEL_EXTERIOR: "#c3c2b7",   # neutral grey (not a category hue)
}
ELEMENT_TEXT = {label: "#ffffff" for label in ELEMENT_COLORS}
ELEMENT_TEXT[LABEL_INTERNAL] = INK          # dark text on yellow
ELEMENT_TEXT[LABEL_EXTERIOR] = INK          # dark text on light grey
ELEMENT_SYMBOLS = {
    LABEL_STEM: "circle",
    LABEL_HAIRPIN: "square",
    LABEL_BULGE: "diamond",
    LABEL_INTERNAL: "triangle-up",
    LABEL_MULTI: "pentagon",
    LABEL_EXTERIOR: "circle",
}

PAIR_LINE = INK_SECONDARY
RETAINED = INK_MUTED
LOST = "#d03b3b"
NEW = "#0ca30c"
MUTATION_MARK = INK
NEGATIVE_POLE = "#2a78d6"      # stabilising
POSITIVE_POLE = "#e34948"      # destabilising
DIVERGING_SCALE = [[0.0, "#184f95"], [0.25, "#6da7ec"], [0.5, "#f0efec"],
                   [0.75, "#ec8a89"], [1.0, "#b52f2f"]]
SEQUENTIAL_SCALE = [[0.0, "#f0efec"], [0.15, "#cde2fb"], [0.5, "#6da7ec"],
                    [0.8, "#256abf"], [1.0, "#0d366b"]]


# ---------------------------- helpers ---------------------------------- #
def _marker_size(n: int) -> float:
    """Large enough to hold a letter for short RNAs, smaller for long ones."""
    return max(7.0, min(20.0, 22.0 - n / 8.0))


def _show_letters(n: int) -> bool:
    return n <= 120


def _base_layout(fig: go.Figure, title: str, height: int) -> None:
    fig.update_layout(
        title=dict(text=title, x=0.0, xanchor="left", font=dict(size=16, color=INK)),
        height=height,
        paper_bgcolor=SURFACE,
        plot_bgcolor=SURFACE,
        font=dict(family=FONT, color=INK_SECONDARY, size=13),
        margin=dict(l=20, r=20, t=70, b=40),
        legend=dict(orientation="h", x=0.0, xanchor="left", y=1.0, yanchor="bottom",
                    font=dict(size=12, color=INK_SECONDARY)),
        hoverlabel=dict(bgcolor="#ffffff", bordercolor=GRID, font=dict(family=FONT, color=INK)),
        hovermode="closest",
    )


def _hide_axes(fig: go.Figure, equal_aspect: bool = False) -> None:
    fig.update_xaxes(visible=False, showgrid=False, zeroline=False)
    if equal_aspect:
        fig.update_yaxes(visible=False, showgrid=False, zeroline=False,
                         scaleanchor="x", scaleratio=1)
    else:
        fig.update_yaxes(visible=False, showgrid=False, zeroline=False)


def nucleotide_hover(sequence: str, pairs: List[Pair], elements=None,
                     extra: Optional[Dict[int, str]] = None) -> List[str]:
    """One hover text per nucleotide: position, base, partner, element."""
    n = len(sequence)
    partner = pair_partner_map(pairs, n)
    elements = elements or analyse_structure(sequence, pairs)
    texts = []
    for k in range(n):
        lines = [f"<b>{sequence[k]}{k + 1}</b>"]
        if partner[k] >= 0:
            p = partner[k]
            b1, b2 = (sequence[k], sequence[p]) if k < p else (sequence[p], sequence[k])
            lines.append(f"Paired with {sequence[p]}{p + 1} ({pair_type(b1, b2)})")
        else:
            lines.append("Unpaired")
        lines.append(f"Element: {elements.position_labels[k]} ({elements.position_element[k]})")
        if extra and k in extra:
            lines.append(extra[k])
        texts.append("<br>".join(lines))
    return texts


def _nucleotide_traces(fig: go.Figure, sequence: str, xs: Sequence[float], ys: Sequence[float],
                       elements, hover: List[str]) -> None:
    """One marker trace per element type, so the legend names each element."""
    n = len(sequence)
    size = _marker_size(n)
    letters = _show_letters(n)
    for label in ELEMENT_ORDER:
        idx = [k for k in range(n) if elements.position_labels[k] == label]
        if not idx:
            continue
        fig.add_trace(go.Scatter(
            x=[xs[k] for k in idx], y=[ys[k] for k in idx],
            mode="markers+text" if letters else "markers",
            text=[sequence[k] for k in idx] if letters else None,
            textposition="middle center",
            textfont=dict(color=ELEMENT_TEXT[label], size=max(8, size * 0.55), family=FONT),
            marker=dict(size=size, color=ELEMENT_COLORS[label], symbol=ELEMENT_SYMBOLS[label],
                        line=dict(color=SURFACE, width=1.5)),
            hovertext=[hover[k] for k in idx], hoverinfo="text",
            name=label, legendgroup=label,
        ))


def _line_trace(segments: List[Tuple[List[float], List[float]]], color: str, width: float,
                name: str, dash: str = "solid", show_legend: bool = True,
                opacity: float = 1.0) -> go.Scatter:
    """Many line segments in ONE trace (separated by None) - fast for long RNAs."""
    xs: List[Optional[float]] = []
    ys: List[Optional[float]] = []
    for sx, sy in segments:
        xs.extend(list(sx) + [None])
        ys.extend(list(sy) + [None])
    return go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=color, width=width, dash=dash),
                      hoverinfo="skip", name=name, showlegend=show_legend and bool(segments),
                      opacity=opacity)


def _split_by_pair_type(sequence: str, pairs: List[Pair]) -> Tuple[List[Pair], List[Pair]]:
    canonical = [p for p in pairs if pair_type(sequence[p[0]], sequence[p[1]]) == "Watson-Crick"]
    wobble = [p for p in pairs if pair_type(sequence[p[0]], sequence[p[1]]) == "Wobble"]
    return canonical, wobble


# ------------------------------------------------------------------ #
# 1. Arc diagram
# ------------------------------------------------------------------ #
def arc_diagram(sequence: str, pairs: List[Pair], title: str = "Arc diagram",
                height: int = 420) -> go.Figure:
    n = len(sequence)
    elements = analyse_structure(sequence, pairs)
    fig = go.Figure()
    canonical, wobble = _split_by_pair_type(sequence, pairs)
    fig.add_trace(_line_trace([arc_points(i, j) for i, j in canonical], PAIR_LINE, 1.6,
                              "Watson-Crick pair"))
    fig.add_trace(_line_trace([arc_points(i, j) for i, j in wobble], PAIR_LINE, 1.6,
                              "G-U wobble pair", dash="dash"))
    xs = list(range(1, n + 1))
    _nucleotide_traces(fig, sequence, xs, [0.0] * n, elements, nucleotide_hover(sequence, pairs, elements))

    max_radius = max([(j - i) / 2 for i, j in pairs], default=1.0)
    _base_layout(fig, title, height)
    fig.update_xaxes(title=dict(text="Nucleotide position", font=dict(color=INK_MUTED)),
                     range=[0, n + 1], showgrid=False, zeroline=False,
                     showline=True, linecolor=BASELINE, tickfont=dict(color=INK_MUTED))
    fig.update_yaxes(visible=False, range=[-max(1.5, max_radius * 0.08), max_radius * 1.08 + 1])
    return fig


# ------------------------------------------------------------------ #
# 2. Circle plot
# ------------------------------------------------------------------ #
def circle_plot(sequence: str, pairs: List[Pair], title: str = "Circle plot",
                height: int = 520) -> go.Figure:
    n = len(sequence)
    elements = analyse_structure(sequence, pairs)
    coords = circle_positions(n)
    fig = go.Figure()
    # backbone
    fig.add_trace(go.Scatter(x=[c[0] for c in coords], y=[c[1] for c in coords], mode="lines",
                             line=dict(color=BASELINE, width=1), hoverinfo="skip",
                             name="Backbone", showlegend=False))
    canonical, wobble = _split_by_pair_type(sequence, pairs)
    fig.add_trace(_line_trace([chord_points(coords[i], coords[j]) for i, j in canonical],
                              PAIR_LINE, 1.4, "Watson-Crick pair"))
    fig.add_trace(_line_trace([chord_points(coords[i], coords[j]) for i, j in wobble],
                              PAIR_LINE, 1.4, "G-U wobble pair", dash="dash"))
    _nucleotide_traces(fig, sequence, [c[0] for c in coords], [c[1] for c in coords],
                       elements, nucleotide_hover(sequence, pairs, elements))
    fig.add_annotation(x=coords[0][0] * 1.13, y=coords[0][1] * 1.13, text="5'", showarrow=False,
                       font=dict(color=INK_SECONDARY, size=13))
    fig.add_annotation(x=coords[-1][0] * 1.13, y=coords[-1][1] * 1.13, text="3'", showarrow=False,
                       font=dict(color=INK_SECONDARY, size=13))
    _base_layout(fig, title, height)
    _hide_axes(fig, equal_aspect=True)
    return fig


# ------------------------------------------------------------------ #
# 3. Planar schematic
# ------------------------------------------------------------------ #
def planar_plot(sequence: str, pairs: List[Pair],
                title: str = "Schematic secondary structure (2-D, not a 3-D model)",
                height: int = 560,
                pair_status: Optional[Dict[Pair, str]] = None,
                highlight: Optional[int] = None,
                highlight_label: str = "Mutated nucleotide",
                coords: Optional[List[Tuple[float, float]]] = None,
                extra_hover: Optional[Dict[int, str]] = None,
                lost_pairs: Optional[List[Pair]] = None) -> go.Figure:
    """
    pair_status: optional {pair: "retained" | "new"} used in mutant comparisons.
    lost_pairs:  optional pairs drawn red-dashed (wild-type pairs lost in a mutant).
    highlight:   optional 0-based position drawn with a ring (the mutation).
    """
    n = len(sequence)
    elements = analyse_structure(sequence, pairs)
    coords = coords or planar_layout(n, pairs)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[c[0] for c in coords], y=[c[1] for c in coords], mode="lines",
                             line=dict(color=BASELINE, width=1.5), hoverinfo="skip",
                             name="Backbone", showlegend=False))

    def segment(p: Pair):
        return [coords[p[0]][0], coords[p[1]][0]], [coords[p[0]][1], coords[p[1]][1]]

    if pair_status is None:
        canonical, wobble = _split_by_pair_type(sequence, pairs)
        fig.add_trace(_line_trace([segment(p) for p in canonical], PAIR_LINE, 2.0, "Watson-Crick pair"))
        fig.add_trace(_line_trace([segment(p) for p in wobble], PAIR_LINE, 2.0, "G-U wobble pair",
                                  dash="dash"))
    else:
        retained = [p for p in pairs if pair_status.get(p) == "retained"]
        new = [p for p in pairs if pair_status.get(p) == "new"]
        fig.add_trace(_line_trace([segment(p) for p in retained], RETAINED, 2.0, "Retained pair"))
        fig.add_trace(_line_trace([segment(p) for p in new], NEW, 3.0, "New pair in mutant"))
    if lost_pairs:
        fig.add_trace(_line_trace([segment(p) for p in lost_pairs], LOST, 3.0, "Lost in mutant",
                                  dash="dash"))

    hover = nucleotide_hover(sequence, pairs, elements, extra_hover)
    _nucleotide_traces(fig, sequence, [c[0] for c in coords], [c[1] for c in coords], elements, hover)

    if highlight is not None:
        fig.add_trace(go.Scatter(
            x=[coords[highlight][0]], y=[coords[highlight][1]], mode="markers",
            marker=dict(size=_marker_size(n) + 12, color="rgba(0,0,0,0)", symbol="circle",
                        line=dict(color=MUTATION_MARK, width=3)),
            hovertext=[hover[highlight]], hoverinfo="text", name=highlight_label,
        ))
    fig.add_annotation(x=coords[0][0], y=coords[0][1], text="5'", showarrow=False,
                       xshift=-18, font=dict(color=INK_SECONDARY, size=13))
    _base_layout(fig, title, height)
    _hide_axes(fig, equal_aspect=True)
    return fig


# ------------------------------------------------------------------ #
# 4. Mutation comparison
# ------------------------------------------------------------------ #
def mutation_arc_comparison(analysis, height: int = 520) -> go.Figure:
    """Wild-type pairs above the axis, mutant pairs below; changes highlighted."""
    wt = analysis.wild_type
    mut = analysis.mutant
    comp = analysis.comparison
    m = analysis.mutation
    n = len(wt.sequence)
    fig = go.Figure()

    fig.add_trace(_line_trace([arc_points(i, j, +1) for i, j in comp.retained], RETAINED, 1.6,
                              "Retained pair"))
    fig.add_trace(_line_trace([arc_points(i, j, -1) for i, j in comp.retained], RETAINED, 1.6,
                              "Retained pair", show_legend=False))
    fig.add_trace(_line_trace([arc_points(i, j, +1) for i, j in comp.lost], LOST, 2.4,
                              "Lost in mutant", dash="dash"))
    fig.add_trace(_line_trace([arc_points(i, j, -1) for i, j in comp.new], NEW, 2.4,
                              "New in mutant"))

    wt_partner = pair_partner_map(wt.pairs, n)
    mut_partner = pair_partner_map(mut.pairs, n)
    hover = []
    for k in range(n):
        wt_text = f"paired with {wt_partner[k] + 1}" if wt_partner[k] >= 0 else "unpaired"
        mut_text = f"paired with {mut_partner[k] + 1}" if mut_partner[k] >= 0 else "unpaired"
        base = wt.sequence[k] if k != m.index else f"{m.original}→{m.replacement}"
        hover.append(f"<b>Position {k + 1}: {base}</b><br>Wild type: {wt_text}<br>Mutant: {mut_text}")

    size = _marker_size(n)
    letters = _show_letters(n)
    changed = set(comp.changed_positions)
    colors = [INK if k == m.index else (INK_SECONDARY if k in changed else BASELINE) for k in range(n)]
    fig.add_trace(go.Scatter(
        x=list(range(1, n + 1)), y=[0.0] * n,
        mode="markers+text" if letters else "markers",
        text=[mut.sequence[k] for k in range(n)] if letters else None,
        textposition="middle center",
        textfont=dict(color=[("#ffffff" if c != BASELINE else INK) for c in colors],
                      size=max(8, size * 0.55)),
        marker=dict(size=size, color=colors, line=dict(color=SURFACE, width=1.5)),
        hovertext=hover, hoverinfo="text", name="Nucleotide (dark = pairing changed)",
    ))
    fig.add_trace(go.Scatter(
        x=[m.position], y=[0.0], mode="markers",
        marker=dict(size=size + 12, color="rgba(0,0,0,0)", line=dict(color=MUTATION_MARK, width=3)),
        hovertext=[hover[m.index]], hoverinfo="text", name=f"Mutation {m.label}",
    ))

    radius = max([(j - i) / 2 for i, j in wt.pairs + mut.pairs], default=1.0)
    fig.add_annotation(x=0.5, y=radius * 0.95, text="Wild type", showarrow=False, xanchor="left",
                       font=dict(color=INK, size=13))
    fig.add_annotation(x=0.5, y=-radius * 0.95, text="Mutant", showarrow=False, xanchor="left",
                       font=dict(color=INK, size=13))
    _base_layout(fig, f"Wild type vs mutant {m.label} ({analysis.algorithm})", height)
    fig.update_xaxes(title=dict(text="Nucleotide position", font=dict(color=INK_MUTED)),
                     range=[0, n + 1], showgrid=False, zeroline=False, tickfont=dict(color=INK_MUTED))
    fig.update_yaxes(visible=False, range=[-radius * 1.1 - 1, radius * 1.1 + 1])
    fig.add_hline(y=0, line_color=BASELINE, line_width=1)
    return fig


def mutation_planar_plots(analysis, height: int = 520) -> Tuple[go.Figure, go.Figure]:
    """Two schematic drawings (wild type, mutant) with the mutated base ringed."""
    wt, mut, comp, m = analysis.wild_type, analysis.mutant, analysis.comparison, analysis.mutation
    wt_status = {p: "retained" for p in comp.retained}
    mut_status = {p: "retained" for p in comp.retained}
    mut_status.update({p: "new" for p in comp.new})

    # Wild type: retained pairs grey, lost pairs red dashed.
    fig_wt = planar_plot(wt.sequence, wt.pairs, title=f"Wild type ({m.original}{m.position})",
                         height=height, pair_status=wt_status, lost_pairs=comp.lost,
                         highlight=m.index, highlight_label="Mutation site")
    fig_mut = planar_plot(mut.sequence, mut.pairs, title=f"Mutant ({m.label})", height=height,
                          pair_status=mut_status, highlight=m.index, highlight_label="Mutated nucleotide")
    return fig_wt, fig_mut


def mutation_scan_heatmap(rows: List[dict], sequence: str, height: int = 330) -> go.Figure:
    """Rows from mutation_analysis.scan_mutations(); ddG by position (x) and new base (y)."""
    n = len(sequence)
    bases = ["A", "U", "G", "C"]
    z = [[None] * n for _ in bases]
    text = [[f"{sequence[k]}{k + 1}: wild-type base" for k in range(n)] for _ in bases]
    for r in rows:
        y = bases.index(r["Replacement"])
        x = r["Position"] - 1
        z[y][x] = r["ΔΔG (kcal/mol)"]
        text[y][x] = (f"<b>{r['Mutation']}</b><br>ΔΔG = {r['ΔΔG (kcal/mol)']:+.2f} kcal/mol"
                      f"<br>Pairs lost {r['Pairs lost']}, gained {r['Pairs gained']}")
    limit = max([abs(r["ΔΔG (kcal/mol)"]) for r in rows] + [0.5])
    fig = go.Figure(go.Heatmap(
        z=z, x=list(range(1, n + 1)), y=[f"→ {b}" for b in bases],
        colorscale=DIVERGING_SCALE, zmid=0, zmin=-limit, zmax=limit,
        text=text, hovertemplate="%{text}<extra></extra>", xgap=1, ygap=1,
        colorbar=dict(title=dict(text="ΔΔG<br>(kcal/mol)"), thickness=12),
    ))
    _base_layout(fig, "Mutation scan: ΔΔG of every single-nucleotide mutant", height)
    fig.update_xaxes(title=dict(text="Position (wild-type base shown in hover; blank = no change)"),
                     showgrid=False, tickfont=dict(color=INK_MUTED))
    fig.update_yaxes(title=dict(text="New base"), showgrid=False, tickfont=dict(color=INK_MUTED))
    return fig


# ------------------------------------------------------------------ #
# 5. Energy and DP teaching figures
# ------------------------------------------------------------------ #
def loop_energy_bar(sequence: str, pairs: List[Pair], height: Optional[int] = None) -> go.Figure:
    """Horizontal bar per loop: negative = stabilising, positive = destabilising."""
    evaluation = evaluate_structure(sequence, pairs)
    loops = [l for l in evaluation.loops if not (l.loop_type == "Exterior" and l.energy == 0)]
    labels = [f"{l.describe(sequence)}" for l in loops]
    values = [round(l.energy, 2) for l in loops]
    fig = go.Figure()
    for name, color, keep in (("Stabilising (ΔG < 0)", NEGATIVE_POLE, lambda v: v < 0),
                              ("Destabilising (ΔG ≥ 0)", POSITIVE_POLE, lambda v: v >= 0)):
        idx = [k for k, v in enumerate(values) if keep(v)]
        if idx:
            fig.add_trace(go.Bar(
                x=[values[k] for k in idx], y=[labels[k] for k in idx], orientation="h",
                marker=dict(color=color, line=dict(color=SURFACE, width=1)), name=name,
                hovertemplate="%{y}<br>%{x:+.2f} kcal/mol<extra></extra>",
            ))
    total = round(evaluation.total, 2)
    _base_layout(fig, f"Energy contribution of each loop (total {total:.2f} kcal/mol)",
                 height or max(260, 28 * len(loops) + 120))
    fig.update_layout(barmode="relative", bargap=0.25)
    fig.update_xaxes(title=dict(text="Estimated ΔG contribution (kcal/mol)"), showgrid=True,
                     gridcolor=GRID, zeroline=True, zerolinecolor=BASELINE)
    fig.update_yaxes(autorange="reversed", showgrid=False, categoryorder="array", categoryarray=labels)
    return fig


def dp_matrix_heatmap(matrix, sequence: str, title: str = "Nussinov DP matrix N(i, j)",
                      height: int = 520) -> go.Figure:
    """Upper triangle of the DP table; each cell = max pairs in sequence[i..j]."""
    n = len(sequence)
    z = [[(float(matrix[i][j]) if j >= i else None) for j in range(n)] for i in range(n)]
    labels = [f"{sequence[k]}{k + 1}" for k in range(n)]
    fig = go.Figure(go.Heatmap(
        z=z, x=labels, y=labels, colorscale=SEQUENTIAL_SCALE,
        hovertemplate="i = %{y}, j = %{x}<br>N(i, j) = %{z}<extra></extra>",
        colorbar=dict(title=dict(text="Pairs"), thickness=12),
    ))
    _base_layout(fig, title, height)
    fig.update_xaxes(title=dict(text="j (end position)"), showgrid=False, side="top",
                     tickfont=dict(color=INK_MUTED, size=10))
    fig.update_yaxes(title=dict(text="i (start position)"), autorange="reversed", showgrid=False,
                     tickfont=dict(color=INK_MUTED, size=10), scaleanchor="x", scaleratio=1)
    return fig
