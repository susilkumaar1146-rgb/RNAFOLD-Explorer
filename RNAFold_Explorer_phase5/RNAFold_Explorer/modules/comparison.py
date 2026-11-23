"""
comparison.py
=============
Side-by-side comparison of the Nussinov and thermodynamic predictions for
the same sequence.

Why can the two algorithms disagree?
------------------------------------
* Nussinov treats every base pair as equally good (+1) and ignores loops,
  so it happily adds isolated pairs and long, costly loops.
* The thermodynamic model scores LOOPS: stacked pairs are rewarded, loops
  are penalised. It may leave bases unpaired if pairing them would cost
  more loop energy than the stacking it gains.
To compare fairly, the Nussinov structure is also scored with the same
energy model (energy_eval) - this is a separate evaluation, not part of the
Nussinov algorithm itself.
"""

from dataclasses import dataclass
from typing import List, Tuple

from modules.energy_eval import evaluate_structure
from modules.fold_result import FoldResult
from modules.structure_elements import analyse_structure

Pair = Tuple[int, int]


@dataclass
class AlgorithmComparison:
    nussinov: FoldResult
    thermodynamic: FoldResult
    nussinov_energy: float           # Nussinov structure scored with the energy model
    shared: List[Pair]
    only_nussinov: List[Pair]
    only_thermodynamic: List[Pair]
    position_agreement: float        # % of positions with the same partner (or both unpaired)

    @property
    def identical(self) -> bool:
        return not self.only_nussinov and not self.only_thermodynamic

    def table_rows(self) -> List[dict]:
        n_el = analyse_structure(self.nussinov.sequence, self.nussinov.pairs).counts()
        t_el = analyse_structure(self.thermodynamic.sequence, self.thermodynamic.pairs).counts()
        nus, th = self.nussinov, self.thermodynamic
        rows = [
            ("Prediction objective", "Maximum number of base pairs", "Minimum estimated free energy"),
            ("Number of base pairs", str(nus.pair_count), str(th.pair_count)),
            ("Dot-bracket notation", nus.dot_bracket, th.dot_bracket),
            ("Estimated free energy", "Not calculated by Nussinov", f"{th.energy:.2f} kcal/mol"),
            ("Energy of this structure under the thermodynamic model (separate evaluation)",
             f"{self.nussinov_energy:.2f} kcal/mol", f"{th.energy:.2f} kcal/mol"),
            ("Stacking energy", "Not considered", "Considered (nearest-neighbour stacks)"),
            ("Loop penalties", "Only the minimum hairpin size (3 nt)",
             "Hairpin, bulge, internal and multiloop penalties"),
            ("Equally optimal structures",
             (f"{nus.optimal_structure_count:,}" if nus.optimal_structure_count else "not counted"),
             "Not enumerated (one minimum shown)"),
        ]
        for key in ("Stems", "Hairpin loops", "Bulge loops", "Internal loops", "Multiloops"):
            rows.append((key, str(n_el[key]), str(t_el[key])))
        return [{"Parameter": p, "Nussinov": a, "Thermodynamic model": b} for p, a, b in rows]

    def explanation(self) -> List[str]:
        nus, th = self.nussinov, self.thermodynamic
        lines = []
        if self.identical:
            lines.append("Both algorithms predict exactly the same base pairs for this sequence.")
        else:
            lines.append(
                f"The predictions share {len(self.shared)} base pair(s); Nussinov has "
                f"{len(self.only_nussinov)} pair(s) not in the thermodynamic structure, and the "
                f"thermodynamic structure has {len(self.only_thermodynamic)} pair(s) not in the "
                f"Nussinov structure. {self.position_agreement:.0f}% of positions have the same "
                "pairing status in both.")
        if nus.pair_count > th.pair_count:
            lines.append(
                f"Nussinov forms more pairs ({nus.pair_count} vs {th.pair_count}) because every pair "
                "counts +1. Under the energy model its structure scores "
                f"{self.nussinov_energy:.2f} kcal/mol versus {th.energy:.2f} kcal/mol: the extra pairs "
                "create loops whose penalties outweigh their stacking.")
        elif self.nussinov_energy > th.energy and not self.identical:
            lines.append(
                f"Under the energy model the Nussinov structure scores {self.nussinov_energy:.2f} "
                f"kcal/mol, less favourable than the thermodynamic prediction "
                f"({th.energy:.2f} kcal/mol), because stacking and loop costs differ.")
        lines.append(
            "By construction, the thermodynamic prediction has the lowest estimated energy of all "
            "structures whose bulges and internal loops have at most 30 unpaired nucleotides, "
            "while the Nussinov prediction has the largest possible number of pairs.")
        return lines


def compare_predictions(nussinov_result: FoldResult, thermo_result: FoldResult) -> AlgorithmComparison:
    if nussinov_result.sequence != thermo_result.sequence:
        raise ValueError("Both predictions must be for the same sequence.")
    seq = nussinov_result.sequence
    a, b = set(nussinov_result.pairs), set(thermo_result.pairs)
    partner_a, partner_b = [-1] * len(seq), [-1] * len(seq)
    for i, j in a:
        partner_a[i], partner_a[j] = j, i
    for i, j in b:
        partner_b[i], partner_b[j] = j, i
    same = sum(1 for k in range(len(seq)) if partner_a[k] == partner_b[k])
    return AlgorithmComparison(
        nussinov=nussinov_result,
        thermodynamic=thermo_result,
        nussinov_energy=round(evaluate_structure(seq, nussinov_result.pairs).total, 2) + 0.0,
        shared=sorted(a & b),
        only_nussinov=sorted(a - b),
        only_thermodynamic=sorted(b - a),
        position_agreement=100.0 * same / len(seq) if seq else 100.0,
    )
