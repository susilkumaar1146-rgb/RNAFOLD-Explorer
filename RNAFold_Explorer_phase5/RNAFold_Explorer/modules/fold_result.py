"""
fold_result.py
==============
A small container that every prediction algorithm returns, so the rest of
the application (tables, plots, mutation analysis, downloads) can treat
Nussinov and thermodynamic results in the same way.

Position convention
-------------------
Inside the code, positions are 0-based (Python indexing).
Anything shown to the user is 1-based (biological numbering).
pairs_1based() performs that conversion in a single place.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from modules.pairing import pair_type


@dataclass
class FoldResult:
    sequence: str
    algorithm: str                           # "Nussinov" or "Thermodynamic"
    pairs: List[Tuple[int, int]]             # 0-based (i, j) with i < j, sorted
    dot_bracket: str
    pair_count: int
    energy: Optional[float] = None           # kcal/mol; None for Nussinov
    alternative_optima: bool = False         # other structures with equal score exist
    optimal_structure_count: Optional[int] = None   # Nussinov: number of co-optimal structures
    notes: List[str] = field(default_factory=list)

    def pairs_1based(self) -> List[Tuple[int, int]]:
        return [(i + 1, j + 1) for i, j in self.pairs]

    def pair_table(self) -> List[dict]:
        """Rows describing every base pair, ready for a pandas DataFrame."""
        rows = []
        for i, j in self.pairs:
            b1, b2 = self.sequence[i], self.sequence[j]
            rows.append({
                "Position i": i + 1,
                "Base i": b1,
                "Position j": j + 1,
                "Base j": b2,
                "Pair": f"{b1}-{b2}",
                "Type": pair_type(b1, b2),
            })
        return rows
