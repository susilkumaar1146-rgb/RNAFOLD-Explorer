"""
thermodynamic.py
================
Simplified thermodynamic RNA folding: a minimum-free-energy (MFE) dynamic
programming algorithm in the style of Zuker & Stiegler (1981), using the
simplified Turner-2004-based energy model in energy_parameters.py.

This is NOT a full Zuker/mfold or ViennaRNA implementation. Terminal
mismatches, dangling ends, special hairpins, small internal-loop tables
and coaxial stacking are not modelled (see data/energy_parameters.json).

Difference from Nussinov
------------------------
Nussinov gives every pair a score of 1. Here a pair by itself has NO
energy: energy comes from LOOPS. Stacked pairs are stabilising (negative),
loops cost energy (positive). Stacking depends on TWO neighbouring pairs,
so the DP must remember which pair closes each region. That is why we
need the extra state V(i, j).

Dynamic programming states (0-based indices)
--------------------------------------------
V[i, j]  : minimum energy of the region i..j GIVEN THAT i pairs with j.
           The pair (i, j) closes exactly one loop, which can be:
             (a) a hairpin:            H(i, j)
             (b) a stack/bulge/internal loop with ONE inner pair (p, q):
                                       interior(i, j, p, q) + V[p, q]
                 (at most MAX_LOOP = 30 unpaired bases in the loop)
             (c) a multiloop with two or more inner helices:
                                       a + c + AU/GU(i,j)
                                       + min_k WM[i+1, k] + WM[k+1, j-1]

WM[i, j] : minimum energy of region i..j when it lies INSIDE a multiloop
           and contains at least one helix:
             WM[i, j] = min( V[i, j] + c + AU/GU(i, j)   one helix i..j
                             WM[i+1, j] + b              i unpaired
                             WM[i, j-1] + b              j unpaired
                             min_k WM[i, k] + WM[k+1, j] two parts )

W[j]     : minimum energy of the first j bases (the exterior loop):
             W[0] = 0
             W[j] = min( W[j-1]                          base j-1 unpaired
                         min_i W[i] + V[i, j-1] + AU/GU(i, j-1) )
           The predicted MFE is W[n].

The exterior W is 1-dimensional because the exterior loop always starts
at the 5' end. This is the same as W(1, j) in the 2-D notation.

Traceback
---------
For each cell we store a back-pointer saying which option gave the minimum.
Following the pointers from W[n] recovers the pairs.

Complexity: O(n^3) for the multiloop/bifurcation terms plus
O(n^2 * MAX_LOOP^2) for bulges and internal loops; O(n^2) memory.
"""

import math
from typing import Dict, List, Optional, Tuple

import numpy as np

from modules.dot_bracket import pairs_to_dot_bracket
from modules.energy_parameters import EnergyModel
from modules.fold_result import FoldResult
from modules.pairing import MIN_HAIRPIN, can_pair

INF = math.inf
# Measured: ~1 s at 200 nt, ~3 s at 300 nt (pure Python; a laptop may be slower)
MAX_LENGTH = 400       # hard limit for the user interface
WARN_LENGTH = 200      # show a "this may take a while" warning above this

_default_model: Optional[EnergyModel] = None


def get_default_model() -> EnergyModel:
    """Load the parameter file once and reuse it."""
    global _default_model
    if _default_model is None:
        _default_model = EnergyModel()
    return _default_model


class ThermoTables:
    """Holds the DP matrices and back-pointers after filling."""

    def __init__(self, n: int):
        self.n = n
        self.V = np.full((n, n), INF)
        self.WM = np.full((n, n), INF)
        self.W = np.zeros(n + 1)
        self.V_bt: Dict[Tuple[int, int], tuple] = {}
        self.WM_bt: Dict[Tuple[int, int], tuple] = {}
        self.W_bt: List[tuple] = [("start",)] * (n + 1)


def fill_tables(seq: str, model: EnergyModel) -> ThermoTables:
    n = len(seq)
    t = ThermoTables(n)
    V, WM = t.V, t.WM
    max_loop = model.max_loop
    b = model.multiloop_unpaired()

    for span in range(MIN_HAIRPIN + 1, n):            # span = j - i
        for i in range(0, n - span):
            j = i + span

            # ---------------- V(i, j): i pairs with j ---------------- #
            if can_pair(seq[i], seq[j]):
                best = model.hairpin(seq, i, j)        # (a) hairpin
                best_bt: tuple = ("H",)

                # (b) stack / bulge / internal loop with inner pair (p, q)
                p_max = min(i + 1 + max_loop, j - MIN_HAIRPIN - 2)
                for p in range(i + 1, p_max + 1):
                    left = p - i - 1
                    q_min = max(p + MIN_HAIRPIN + 1, j - 1 - (max_loop - left))
                    V_row = V[p]
                    for q in range(j - 1, q_min - 1, -1):
                        inner = V_row[q]
                        if inner == INF:
                            continue
                        e = inner + model.interior(seq, i, j, p, q)
                        if e < best:
                            best, best_bt = e, ("I", p, q)

                # (c) multiloop: split i+1..j-1 into two WM parts at k
                if j - i - 1 >= 2 * (MIN_HAIRPIN + 2):
                    left_parts = WM[i + 1, i + 1:j - 1]     # WM[i+1, k], k = i+1..j-2
                    right_parts = WM[i + 2:j, j - 1]        # WM[k+1, j-1]
                    sums = left_parts + right_parts
                    idx = int(np.argmin(sums))
                    if sums[idx] < INF:
                        e = sums[idx] + model.multiloop_closing(seq, i, j)
                        if e < best:
                            best, best_bt = e, ("M", i + 1 + idx)

                V[i, j] = best
                t.V_bt[(i, j)] = best_bt

            # ------------- WM(i, j): part of a multiloop ------------- #
            best = INF
            best_bt = ("none",)
            if V[i, j] < INF:
                e = V[i, j] + model.multiloop_branch(seq, i, j)
                if e < best:
                    best, best_bt = e, ("P",)
            if WM[i + 1, j] + b < best:
                best, best_bt = WM[i + 1, j] + b, ("L",)
            if WM[i, j - 1] + b < best:
                best, best_bt = WM[i, j - 1] + b, ("R",)
            sums = WM[i, i:j] + WM[i + 1:j + 1, j]        # k = i..j-1
            idx = int(np.argmin(sums))
            if sums[idx] < best:
                best, best_bt = sums[idx], ("B", i + idx)
            WM[i, j] = best
            if best < INF:
                t.WM_bt[(i, j)] = best_bt

    # ---------------- W: the exterior loop ---------------- #
    W = t.W
    for j in range(1, n + 1):
        best, best_bt = W[j - 1], ("U",)                  # base j-1 unpaired
        last = j - 1
        for i in range(0, last - MIN_HAIRPIN):
            if V[i, last] < INF:
                e = W[i] + V[i, last] + model.exterior_branch(seq, i, last)
                if e < best:
                    best, best_bt = e, ("P", i)
        W[j] = best
        t.W_bt[j] = best_bt
    return t


def traceback(t: ThermoTables) -> List[Tuple[int, int]]:
    """Follow back-pointers from W[n] and return the sorted 0-based pairs."""
    pairs: List[Tuple[int, int]] = []
    stack: List[tuple] = []

    j = t.n
    while j > 0:                                   # walk the exterior loop
        bt = t.W_bt[j]
        if bt[0] == "U":
            j -= 1
        else:
            i = bt[1]
            stack.append(("V", i, j - 1))
            j = i

    while stack:
        kind, i, j = stack.pop()
        if kind == "V":
            pairs.append((i, j))
            bt = t.V_bt[(i, j)]
            if bt[0] == "I":
                stack.append(("V", bt[1], bt[2]))
            elif bt[0] == "M":
                k = bt[1]
                stack.append(("WM", i + 1, k))
                stack.append(("WM", k + 1, j - 1))
            # "H": hairpin, nothing inside
        else:  # WM
            bt = t.WM_bt[(i, j)]
            if bt[0] == "P":
                stack.append(("V", i, j))
            elif bt[0] == "L":
                stack.append(("WM", i + 1, j))
            elif bt[0] == "R":
                stack.append(("WM", i, j - 1))
            elif bt[0] == "B":
                k = bt[1]
                stack.append(("WM", i, k))
                stack.append(("WM", k + 1, j))
            else:  # pragma: no cover
                raise RuntimeError(f"Invalid multiloop back-pointer at ({i}, {j})")
    return sorted(pairs)


def fold(sequence: str, model: Optional[EnergyModel] = None) -> FoldResult:
    """Predict the minimum-estimated-energy structure of a cleaned RNA sequence."""
    # Imported here to avoid a circular import (energy_eval uses the same model).
    from modules.energy_eval import evaluate_structure

    model = model or get_default_model()
    n = len(sequence)
    if n > MAX_LENGTH:
        raise ValueError(f"Sequence length {n} exceeds the thermodynamic model limit of {MAX_LENGTH} nt.")

    tables = fill_tables(sequence, model)
    pairs = traceback(tables)
    mfe = float(tables.W[n]) if n > 0 else 0.0

    # Independent check: re-score the traced structure loop by loop.
    check = evaluate_structure(sequence, pairs, model)
    if abs(check.total - mfe) > 1e-6:  # pragma: no cover - would indicate a bug
        raise RuntimeError(f"Energy mismatch: DP {mfe:.4f} vs evaluation {check.total:.4f}")

    notes = [
        "Estimated free energy from a simplified nearest-neighbour model "
        "(subset of Turner 2004 parameters, 37 °C). This is a computational estimate, "
        "not an experimentally measured value.",
    ]
    if not pairs:
        notes.append("No structure with negative estimated free energy was found; "
                     "the unfolded chain (0.00 kcal/mol) is predicted.")

    return FoldResult(
        sequence=sequence,
        algorithm="Thermodynamic",
        pairs=pairs,
        dot_bracket=pairs_to_dot_bracket(pairs, n),
        pair_count=len(pairs),
        energy=round(mfe, 2) + 0.0,     # + 0.0 turns -0.0 into 0.0
        notes=notes,
    )
