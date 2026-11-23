"""
energy_eval.py
==============
Computes the estimated free energy of ANY given nested structure using the
same EnergyModel as the thermodynamic DP, and reports the contribution of
every loop.

Why is this useful?
-------------------
1. Correctness check: the energy of the structure returned by the DP
   traceback must equal the DP's minimum energy. If they differ, there is
   a bug in either the recurrences or the traceback.
2. Comparison: it lets us score the Nussinov structure with the SAME
   energy rules ("how stable would the max-pairs structure be?").
3. Mutation analysis and teaching: users see which loops stabilise
   (negative) and which destabilise (positive) the structure.

Loop decomposition
------------------
Each pair (i, j) closes a loop formed by the bases between i and j that are
not inside any inner pair. The loop type is decided by how many pairs
are DIRECTLY inside (i, j) ("children"):
    0 children -> hairpin loop
    1 child    -> stack, bulge or internal loop
    2+ children-> multiloop
The bases not enclosed by any pair form the exterior loop.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from modules.dot_bracket import validate_pairs
from modules.energy_parameters import EnergyModel

Pair = Tuple[int, int]


@dataclass
class LoopContribution:
    loop_type: str                 # Hairpin / Stack / Bulge / Internal loop / Multiloop / Exterior
    closing_pair: Optional[Pair]   # 0-based; None for the exterior loop
    inner_pairs: List[Pair]
    unpaired: int
    energy: float

    def describe(self, sequence: str) -> str:
        if self.closing_pair is None:
            return f"Exterior loop with {len(self.inner_pairs)} helix(es)"
        i, j = self.closing_pair
        text = f"{self.loop_type} closed by {sequence[i]}{i + 1}-{sequence[j]}{j + 1}"
        if self.loop_type in ("Stack", "Bulge", "Internal loop") and self.inner_pairs:
            p, q = self.inner_pairs[0]
            text += f" / {sequence[p]}{p + 1}-{sequence[q]}{q + 1}"
        if self.unpaired and self.loop_type != "Stack":
            text += f" ({self.unpaired} unpaired nt)"
        return text


@dataclass
class EnergyEvaluation:
    total: float
    loops: List[LoopContribution] = field(default_factory=list)


def decompose_loops(pairs: List[Pair], n: int) -> Tuple[List[Pair], dict]:
    """
    Find the direct children of every pair.
    Returns (exterior_children, children_by_pair).
    Uses a stack, scanning the sequence left to right.
    """
    partner = {}
    for i, j in pairs:
        partner[i] = j
        partner[j] = i
    children = {p: [] for p in pairs}
    exterior: List[Pair] = []
    open_stack: List[Pair] = []
    for pos in range(n):
        if pos in partner and partner[pos] > pos:          # opening base
            pair = (pos, partner[pos])
            (children[open_stack[-1]] if open_stack else exterior).append(pair)
            open_stack.append(pair)
        elif pos in partner:                               # closing base
            open_stack.pop()
    return exterior, children


def evaluate_structure(sequence: str, pairs: List[Pair],
                       model: Optional[EnergyModel] = None) -> EnergyEvaluation:
    """Return the total estimated energy and the per-loop contributions."""
    if model is None:
        from modules.thermodynamic import get_default_model
        model = get_default_model()

    n = len(sequence)
    pairs = sorted(pairs)
    errors = validate_pairs(pairs, n, sequence)
    if errors:
        raise ValueError("Cannot evaluate an invalid structure: " + "; ".join(errors))

    exterior, children = decompose_loops(pairs, n)
    loops: List[LoopContribution] = []

    for (i, j) in pairs:
        inner = children[(i, j)]
        unpaired = (j - i - 1) - sum(q - p + 1 for p, q in inner)
        if not inner:
            energy = model.hairpin(sequence, i, j)
            loops.append(LoopContribution("Hairpin", (i, j), [], unpaired, energy))
        elif len(inner) == 1:
            p, q = inner[0]
            energy = model.interior(sequence, i, j, p, q)
            loop_type = model.interior_loop_type(i, j, p, q)
            loops.append(LoopContribution(loop_type, (i, j), [(p, q)], unpaired, energy))
        else:
            energy = model.multiloop_closing(sequence, i, j)
            energy += sum(model.multiloop_branch(sequence, p, q) for p, q in inner)
            energy += unpaired * model.multiloop_unpaired()
            loops.append(LoopContribution("Multiloop", (i, j), list(inner), unpaired, energy))

    exterior_unpaired = n - sum(q - p + 1 for p, q in exterior)
    ext_energy = sum(model.exterior_branch(sequence, p, q) for p, q in exterior)
    loops.append(LoopContribution("Exterior", None, exterior, exterior_unpaired, ext_energy))

    total = sum(loop.energy for loop in loops)
    return EnergyEvaluation(total=total, loops=loops)
