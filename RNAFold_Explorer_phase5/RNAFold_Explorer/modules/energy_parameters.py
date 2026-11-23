"""
energy_parameters.py
====================
Loads data/energy_parameters.json and turns it into an EnergyModel object
with one small function per loop type.

Why a separate module?
----------------------
The thermodynamic dynamic programming (thermodynamic.py) and the
independent energy evaluator (energy_eval.py) must use EXACTLY the same
energy rules. Both call the methods of the same EnergyModel, so a
structure's energy is computed identically everywhere.

Loop types and their energies (kcal/mol)
----------------------------------------
Every nested structure decomposes into loops, each closed by a base pair:

* Stack (i,j) on (i+1,j-1)          -> nearest-neighbour stacking energy
* Hairpin closed by (i,j)            -> initiation(size) + AU/GU end penalty
* Bulge (unpaired bases on ONE side)
      size 1 -> initiation(1) + stacking of the two pairs (helix continues)
      size>1 -> initiation(size) + AU/GU end penalties of both pairs
* Internal loop (unpaired on BOTH sides)
      initiation(size) + 0.6 * |n1 - n2| + 0.7 per AU/GU closing pair
* Multiloop (3 or more helices)      -> a + b*unpaired + c*helices
                                        + AU/GU end penalties
* Exterior loop                      -> AU/GU end penalty of each helix

Negative values stabilise (stacks); positive values destabilise (loops).
"""

import json
import math
import os
from typing import Dict, List, Optional

from modules.pairing import MIN_HAIRPIN, is_au_gu

DEFAULT_PARAMETER_FILE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "energy_parameters.json",
)


class EnergyModel:
    """Simplified nearest-neighbour energy model (subset of Turner 2004)."""

    def __init__(self, parameter_file: Optional[str] = None):
        path = parameter_file or DEFAULT_PARAMETER_FILE
        with open(path, "r", encoding="utf-8") as handle:
            self.raw = json.load(handle)
        p = self.raw

        self.RT = (p["constants"]["gas_constant_kcal_per_mol_K"]
                   * p["constants"]["temperature_K"])

        self.stack_table = self._expand_stacking(p["stacking"]["values"])
        self.terminal_penalty_value = p["terminal_au_gu_penalty"]["value"]

        self.hairpin_table = {int(k): v for k, v in p["hairpin_initiation"]["values"].items()}
        self.bulge_table = {int(k): v for k, v in p["bulge_initiation"]["values"].items()}
        self.internal_table = {int(k): v for k, v in p["internal_initiation"]["values"].items()}
        for k, v in p["internal_initiation"]["small_loops_approximation"]["values"].items():
            self.internal_table[int(k)] = v

        self.internal_asymmetry = p["internal_asymmetry_per_nt"]["value"]
        self.internal_au_gu_closure = p["internal_au_gu_closure"]["value"]

        self.ml_a = p["multiloop"]["a"]
        self.ml_b = p["multiloop"]["b"]
        self.ml_c = p["multiloop"]["c"]

        self.min_hairpin = p["model_settings"]["min_hairpin_size"]
        self.max_loop = p["model_settings"]["max_internal_loop_size"]
        if self.min_hairpin != MIN_HAIRPIN:
            raise ValueError("min_hairpin_size in the JSON file must match pairing.MIN_HAIRPIN")

    # ------------------------------------------------------------------ #
    # Loading helpers
    # ------------------------------------------------------------------ #
    @staticmethod
    def _expand_stacking(values: Dict[str, float]) -> Dict[str, float]:
        """
        Add symmetric equivalents. 'XY/AB' (5'XY3' / 3'AB5') describes the
        same physical stack as 'BA/YX' read from the other strand.
        Raises an error if the file contains two conflicting values.
        """
        table: Dict[str, float] = {}
        for key, value in values.items():
            top, bottom = key.split("/")
            symmetric = bottom[::-1] + "/" + top[::-1]
            for k in (key, symmetric):
                if k in table and abs(table[k] - value) > 1e-9:
                    raise ValueError(f"Conflicting stacking values for {k}")
                table[k] = value
        return table

    def approximate_parameters(self) -> List[str]:
        """Human-readable list of parameters that are NOT Turner 2004 values."""
        p = self.raw
        return [
            "Internal loops of 2 (1x1) and 3 (1x2) nt: " +
            p["internal_initiation"]["small_loops_approximation"]["note"],
            "Multiloop a, b, c: " + p["multiloop"]["note"],
        ]

    # ------------------------------------------------------------------ #
    # Energy terms
    # ------------------------------------------------------------------ #
    def stack(self, seq: str, i: int, j: int, p: int, q: int) -> float:
        """Stacking energy of outer pair (i,j) on the adjacent inner pair (p,q)."""
        key = seq[i] + seq[p] + "/" + seq[j] + seq[q]
        return self.stack_table[key]

    def terminal_penalty(self, seq: str, i: int, j: int) -> float:
        """+0.45 kcal/mol if (i,j) is an A-U or G-U pair ending a helix."""
        return self.terminal_penalty_value if is_au_gu(seq[i], seq[j]) else 0.0

    def _initiation(self, table: Dict[int, float], size: int,
                    ref_size: int, extrapolate) -> float:
        if size in table:
            return table[size]
        largest = max(table)
        if size > largest:
            return extrapolate(size, table[ref_size])
        raise ValueError(f"No initiation parameter for loop size {size}")

    def hairpin_initiation(self, size: int) -> float:
        if size < self.min_hairpin:
            return math.inf
        return self._initiation(
            self.hairpin_table, size, 9,
            lambda n, g9: g9 + 1.75 * self.RT * math.log(n / 9))

    def bulge_initiation(self, size: int) -> float:
        return self._initiation(
            self.bulge_table, size, 6,
            lambda n, g6: g6 + 1.75 * self.RT * math.log(n / 6))

    def internal_initiation(self, size: int) -> float:
        return self._initiation(
            self.internal_table, size, 6,
            lambda n, g6: g6 + 1.08 * math.log(n / 6))

    def hairpin(self, seq: str, i: int, j: int) -> float:
        """Energy of the hairpin loop closed by pair (i, j)."""
        size = j - i - 1
        if size < self.min_hairpin:
            return math.inf
        return self.hairpin_initiation(size) + self.terminal_penalty(seq, i, j)

    def interior(self, seq: str, i: int, j: int, p: int, q: int) -> float:
        """
        Energy of the loop between outer pair (i,j) and inner pair (p,q),
        where i < p < q < j. Covers stacks, bulges and internal loops.
        """
        left = p - i - 1     # unpaired bases on the 5' side
        right = j - q - 1    # unpaired bases on the 3' side

        if left == 0 and right == 0:                       # stacked pair
            return self.stack(seq, i, j, p, q)

        if left == 0 or right == 0:                        # bulge loop
            size = left + right
            if size == 1:
                # helix continues across a 1-nt bulge: add the stack, no end penalty
                return self.bulge_initiation(1) + self.stack(seq, i, j, p, q)
            return (self.bulge_initiation(size)
                    + self.terminal_penalty(seq, i, j)
                    + self.terminal_penalty(seq, p, q))

        size = left + right                                # internal loop
        energy = self.internal_initiation(size)
        energy += self.internal_asymmetry * abs(left - right)
        energy += self.internal_au_gu_closure * (is_au_gu(seq[i], seq[j]) + is_au_gu(seq[p], seq[q]))
        return energy

    # Multiloop terms (linear model: a + b*unpaired + c*helices)
    def multiloop_closing(self, seq: str, i: int, j: int) -> float:
        """Cost of pair (i,j) closing a multiloop: a + c (its own helix) + end penalty."""
        return self.ml_a + self.ml_c + self.terminal_penalty(seq, i, j)

    def multiloop_branch(self, seq: str, k: int, l: int) -> float:
        """Cost of each inner helix (k,l) inside a multiloop: c + end penalty."""
        return self.ml_c + self.terminal_penalty(seq, k, l)

    def multiloop_unpaired(self) -> float:
        """Cost per unpaired nucleotide inside a multiloop: b."""
        return self.ml_b

    def exterior_branch(self, seq: str, k: int, l: int) -> float:
        """A helix in the exterior loop only pays its AU/GU end penalty."""
        return self.terminal_penalty(seq, k, l)

    @staticmethod
    def interior_loop_type(i: int, j: int, p: int, q: int) -> str:
        left, right = p - i - 1, j - q - 1
        if left == 0 and right == 0:
            return "Stack"
        if left == 0 or right == 0:
            return "Bulge"
        return "Internal loop"
