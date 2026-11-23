"""
structure_elements.py
=====================
Identifies the structural elements of a predicted secondary structure.

Biology
-------
A secondary structure is built from:
  * Stems (helices): runs of stacked base pairs (i,j), (i+1,j-1), ...
  * Hairpin loops:  unpaired bases closed by ONE pair, at the tip of a stem
  * Bulge loops:    unpaired bases on ONE side between two stems
  * Internal loops: unpaired bases on BOTH sides between two stems
  * Multiloops:     a loop from which three or more stems radiate
  * Exterior loop:  bases not enclosed by any pair (5'/3' tails, linkers)

Method
------
The loop decomposition from energy_eval.decompose_loops() tells us, for
every pair, which pairs lie directly inside it ("children"):
    0 children  -> hairpin       1 child -> stack, bulge or internal loop
    2+ children -> multiloop     not enclosed -> exterior
Stems are found by following each pair to its child while the child is
(i+1, j-1), i.e. while the two pairs are stacked.

Output
------
analyse_structure() returns a StructureElements object with
  * stems             list of Stem
  * loops             list of Loop (hairpin/bulge/internal/multiloop/exterior)
  * position_labels   one label per nucleotide ("Stem", "Hairpin loop", ...),
                      used to colour the visualisations
  * position_element  one element id per nucleotide ("S1", "H1", "E", ...)
All positions are 0-based internally; the *_table() helpers return 1-based
positions for display.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

from modules.energy_eval import decompose_loops

Pair = Tuple[int, int]

# Label used for each element type (one place, so plots and tables agree).
LABEL_STEM = "Stem"
LABEL_HAIRPIN = "Hairpin loop"
LABEL_BULGE = "Bulge loop"
LABEL_INTERNAL = "Internal loop"
LABEL_MULTI = "Multiloop"
LABEL_EXTERIOR = "Exterior loop"

ELEMENT_ORDER = [LABEL_STEM, LABEL_HAIRPIN, LABEL_BULGE, LABEL_INTERNAL,
                 LABEL_MULTI, LABEL_EXTERIOR]


@dataclass
class Stem:
    element_id: str           # "S1", "S2", ...
    pairs: List[Pair]         # outermost pair first

    @property
    def length(self) -> int:
        return len(self.pairs)

    @property
    def five_prime_range(self) -> Tuple[int, int]:
        return self.pairs[0][0], self.pairs[-1][0]

    @property
    def three_prime_range(self) -> Tuple[int, int]:
        return self.pairs[-1][1], self.pairs[0][1]


@dataclass
class Loop:
    element_id: str           # "H1", "B1", "I1", "M1" or "E"
    label: str                # one of the LABEL_* constants
    closing_pair: object      # Pair, or None for the exterior loop
    inner_pairs: List[Pair]
    unpaired_positions: List[int]

    @property
    def size(self) -> int:
        return len(self.unpaired_positions)


@dataclass
class StructureElements:
    sequence: str
    stems: List[Stem] = field(default_factory=list)
    loops: List[Loop] = field(default_factory=list)
    position_labels: List[str] = field(default_factory=list)
    position_element: List[str] = field(default_factory=list)

    def counts(self) -> Dict[str, int]:
        """Number of each element type (exterior loop excluded)."""
        result = {"Stems": len(self.stems)}
        for label, key in ((LABEL_HAIRPIN, "Hairpin loops"), (LABEL_BULGE, "Bulge loops"),
                           (LABEL_INTERNAL, "Internal loops"), (LABEL_MULTI, "Multiloops")):
            result[key] = sum(1 for loop in self.loops if loop.label == label)
        return result

    def stem_table(self) -> List[dict]:
        rows = []
        for stem in self.stems:
            a, b = stem.five_prime_range
            c, d = stem.three_prime_range
            rows.append({
                "Element": stem.element_id,
                "Base pairs": stem.length,
                "5' strand": f"{a + 1}-{b + 1}",
                "3' strand": f"{c + 1}-{d + 1}",
                "5' sequence": self.sequence[a:b + 1],
                "3' sequence": self.sequence[c:d + 1],
            })
        return rows

    def loop_table(self) -> List[dict]:
        rows = []
        for loop in self.loops:
            if loop.label == LABEL_EXTERIOR and loop.size == 0:
                continue
            closing = "-" if loop.closing_pair is None else \
                f"{loop.closing_pair[0] + 1}-{loop.closing_pair[1] + 1}"
            rows.append({
                "Element": loop.element_id,
                "Type": loop.label,
                "Closing pair": closing,
                "Unpaired nt": loop.size,
                "Unpaired positions": _format_ranges(loop.unpaired_positions),
            })
        return rows


def _format_ranges(positions: List[int]) -> str:
    """[3,4,5,9] (0-based) -> '4-6, 10' (1-based)."""
    if not positions:
        return "-"
    positions = sorted(positions)
    parts, start, prev = [], positions[0], positions[0]
    for p in positions[1:] + [None]:
        if p is not None and p == prev + 1:
            prev = p
            continue
        parts.append(f"{start + 1}" if start == prev else f"{start + 1}-{prev + 1}")
        if p is not None:
            start = prev = p
    return ", ".join(parts)


def _unpaired_inside(i: int, j: int, children: List[Pair]) -> List[int]:
    """Positions between i and j that are not covered by any child pair's span."""
    covered = set()
    for p, q in children:
        covered.update(range(p, q + 1))
    return [k for k in range(i + 1, j) if k not in covered]


def analyse_structure(sequence: str, pairs: List[Pair]) -> StructureElements:
    n = len(sequence)
    pairs = sorted(pairs)
    exterior, children = decompose_loops(pairs, n)
    result = StructureElements(sequence=sequence,
                               position_labels=[LABEL_EXTERIOR] * n,
                               position_element=["E"] * n)

    # ---- stems: start at pairs that are NOT stacked on an outer pair ----
    stacked_inside = {}                         # pair -> its stacked inner pair
    for (i, j), kids in children.items():
        if len(kids) == 1 and kids[0] == (i + 1, j - 1):
            stacked_inside[(i, j)] = kids[0]
    has_stacked_parent = set(stacked_inside.values())
    for start in pairs:                         # pairs sorted 5'->3' = outer first
        if start in has_stacked_parent:
            continue
        stem_pairs = [start]
        while stem_pairs[-1] in stacked_inside:
            stem_pairs.append(stacked_inside[stem_pairs[-1]])
        stem = Stem(f"S{len(result.stems) + 1}", stem_pairs)
        result.stems.append(stem)
        for p, q in stem_pairs:
            result.position_labels[p] = result.position_labels[q] = LABEL_STEM
            result.position_element[p] = result.position_element[q] = stem.element_id

    # ---- loops closed by a pair ----
    counters = {"H": 0, "B": 0, "I": 0, "M": 0}
    for (i, j) in pairs:
        kids = children[(i, j)]
        unpaired = _unpaired_inside(i, j, kids)
        if not kids:
            prefix, label = "H", LABEL_HAIRPIN
        elif len(kids) == 1:
            if not unpaired:
                continue                        # stacked pair: part of a stem
            p, q = kids[0]
            one_sided = (p == i + 1) or (q == j - 1)
            prefix, label = ("B", LABEL_BULGE) if one_sided else ("I", LABEL_INTERNAL)
        else:
            prefix, label = "M", LABEL_MULTI
        counters[prefix] += 1
        loop = Loop(f"{prefix}{counters[prefix]}", label, (i, j), list(kids), unpaired)
        result.loops.append(loop)
        for k in unpaired:
            result.position_labels[k] = label
            result.position_element[k] = loop.element_id

    # ---- exterior loop ----
    covered = set()
    for p, q in exterior:
        covered.update(range(p, q + 1))
    ext_unpaired = [k for k in range(n) if k not in covered]
    result.loops.append(Loop("E", LABEL_EXTERIOR, None, list(exterior), ext_unpaired))
    return result
