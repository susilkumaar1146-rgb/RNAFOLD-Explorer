"""
mutation_analysis.py
====================
Single-nucleotide mutation impact analysis (the project's main extension).

Biology
-------
A point mutation changes one nucleotide. If that base was paired, the pair
may be lost (e.g. G-C -> A-C cannot pair) or weakened (C -> U next to G
gives a G-U wobble). If it was unpaired, new pairs may become possible.
Either way, the whole predicted fold can rearrange, because the dynamic
programming re-optimises the entire sequence.

Workflow
--------
1. create_mutant()      : build the mutant sequence automatically
2. analyse_mutation()   : fold wild type and mutant with the SAME algorithm
                          and parameters, then compare them:
                            - pairs retained / lost / newly formed
                            - dot-bracket differences
                            - structural element counts
                            - Delta-Delta-G = G(mutant) - G(wild type)
3. interpret()          : plain-language summary built ONLY from the computed
                          values, with no claims about biological function
4. scan_mutations()     : (optional) all 3n single mutants -> Delta-Delta-G per
                          position, to see which positions the structure is
                          most sensitive to

Sign convention: Delta-Delta-G > 0 means the mutant has a HIGHER (less
negative) estimated folding free energy, i.e. it is predicted to be less
stable under this model.
"""

import re
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Tuple

from modules import nussinov, thermodynamic
from modules.dot_bracket import pair_partner_map
from modules.energy_eval import evaluate_structure
from modules.fold_result import FoldResult
from modules.pairing import can_pair, pair_type
from modules.structure_elements import analyse_structure

Pair = Tuple[int, int]
BASES = "AUGC"
ALGORITHMS = ("Nussinov", "Thermodynamic")
ZERO_TOLERANCE = 0.005          # |ddG| below this (kcal/mol) is reported as 0
SCAN_MAX_LENGTH = 80            # mutation scan limit: 3n folds (~20 s at 80 nt measured)


# ------------------------------------------------------------------ #
# Creating mutants
# ------------------------------------------------------------------ #
@dataclass(frozen=True)
class Mutation:
    position: int        # 1-based, as shown to the user
    original: str
    replacement: str

    @property
    def label(self) -> str:
        """Standard notation, e.g. 'C12A'."""
        return f"{self.original}{self.position}{self.replacement}"

    @property
    def index(self) -> int:
        return self.position - 1


def create_mutant(sequence: str, position: int, new_base: str) -> Tuple[str, Mutation]:
    """
    Replace the base at `position` (1-based) with `new_base`.
    Raises ValueError for an invalid position, an invalid base, or a
    "mutation" to the same base.
    """
    new_base = new_base.upper().replace("T", "U")
    if not 1 <= position <= len(sequence):
        raise ValueError(f"Position must be between 1 and {len(sequence)} (got {position}).")
    if new_base not in BASES or len(new_base) != 1:
        raise ValueError(f"Replacement base must be one of A, U, G, C (got '{new_base}').")
    original = sequence[position - 1]
    if original == new_base:
        raise ValueError(f"Position {position} is already {original}; choose a different base.")
    mutant = sequence[:position - 1] + new_base + sequence[position:]
    return mutant, Mutation(position, original, new_base)


def parse_mutation(text: str, sequence: Optional[str] = None) -> Mutation:
    """
    Parse standard notation such as 'C12A' (original, position, replacement).
    If a sequence is given, the original base is checked against it.
    """
    match = re.fullmatch(r"\s*([AUGCT])(\d+)([AUGCT])\s*", text.upper())
    if not match:
        raise ValueError("Use the format <original><position><new>, e.g. C12A.")
    original, position, new = match.group(1), int(match.group(2)), match.group(3)
    original, new = original.replace("T", "U"), new.replace("T", "U")
    if sequence is not None:
        if not 1 <= position <= len(sequence):
            raise ValueError(f"Position {position} is outside the sequence (1-{len(sequence)}).")
        if sequence[position - 1] != original:
            raise ValueError(f"Position {position} is {sequence[position - 1]}, not {original}.")
    if original == new:
        raise ValueError("The original and replacement bases are the same.")
    return Mutation(position, original, new)


# ------------------------------------------------------------------ #
# Comparing two structures
# ------------------------------------------------------------------ #
@dataclass
class StructureComparison:
    retained: List[Pair]
    lost: List[Pair]
    new: List[Pair]
    changed_positions: List[int]       # 0-based positions whose dot-bracket symbol OR partner changed

    @property
    def identical(self) -> bool:
        return not self.lost and not self.new


def compare_pairs(wt_pairs: List[Pair], mut_pairs: List[Pair], n: int) -> StructureComparison:
    """Set comparison of two pair lists on sequences of the same length n."""
    wt, mut = set(wt_pairs), set(mut_pairs)
    wt_partner = pair_partner_map(sorted(wt), n)
    mut_partner = pair_partner_map(sorted(mut), n)
    changed = [k for k in range(n) if wt_partner[k] != mut_partner[k]]
    return StructureComparison(
        retained=sorted(wt & mut),
        lost=sorted(wt - mut),
        new=sorted(mut - wt),
        changed_positions=changed,
    )


# ------------------------------------------------------------------ #
# Full analysis
# ------------------------------------------------------------------ #
@dataclass
class MutationAnalysis:
    algorithm: str
    mutation: Mutation
    wild_type: FoldResult
    mutant: FoldResult
    comparison: StructureComparison
    wt_energy: float                 # estimated kcal/mol (see energy_note)
    mut_energy: float
    ddg: float                       # mut_energy - wt_energy
    energy_note: str
    wt_elements: dict = field(default_factory=dict)
    mut_elements: dict = field(default_factory=dict)
    wt_site: str = ""                # description of the mutated position in WT
    mut_site: str = ""

    def summary_rows(self) -> List[dict]:
        """Table of the main quantities, ready for a DataFrame / CSV."""
        rows = [
            ("Mutation", self.mutation.label, ""),
            ("Algorithm", self.algorithm, ""),
            ("Base pairs (wild type)", self.wild_type.pair_count, ""),
            ("Base pairs (mutant)", self.mutant.pair_count, ""),
            ("Pairs retained", len(self.comparison.retained), ""),
            ("Pairs lost", len(self.comparison.lost), ""),
            ("Pairs newly formed", len(self.comparison.new), ""),
            ("Positions with changed pairing", len(self.comparison.changed_positions), ""),
            ("Estimated ΔG wild type", f"{self.wt_energy:.2f}", "kcal/mol"),
            ("Estimated ΔG mutant", f"{self.mut_energy:.2f}", "kcal/mol"),
            ("ΔΔG (mutant − wild type)", f"{self.ddg:+.2f}", "kcal/mol"),
        ]
        return [{"Quantity": q, "Value": v, "Unit": u} for q, v, u in rows]


def _fold(sequence: str, algorithm: str) -> FoldResult:
    if algorithm == "Nussinov":
        return nussinov.fold(sequence)
    if algorithm == "Thermodynamic":
        return thermodynamic.fold(sequence)
    raise ValueError(f"Unknown algorithm '{algorithm}'. Use one of {ALGORITHMS}.")


def _energy(result: FoldResult) -> float:
    """MFE for the thermodynamic model; model energy of the structure for Nussinov."""
    if result.energy is not None:
        return result.energy
    return round(evaluate_structure(result.sequence, result.pairs).total, 2)


def describe_site(result: FoldResult, index: int) -> str:
    """'paired with C9 (G-C, Watson-Crick)' or 'unpaired (Hairpin loop)'."""
    partner = pair_partner_map(result.pairs, len(result.sequence))[index]
    base = result.sequence[index]
    if partner >= 0:
        other = result.sequence[partner]
        b1, b2 = (base, other) if index < partner else (other, base)
        return f"paired with {other}{partner + 1} ({b1}-{b2}, {pair_type(b1, b2)})"
    label = analyse_structure(result.sequence, result.pairs).position_labels[index]
    return f"unpaired ({label.lower()})"


def analyse_mutation(sequence: str, position: int, new_base: str,
                     algorithm: str = "Thermodynamic",
                     wild_type: Optional[FoldResult] = None) -> MutationAnalysis:
    """
    Fold wild type and mutant with the same algorithm and compare them.
    `wild_type` may be passed in to avoid folding the same sequence twice.
    """
    mutant_seq, mutation = create_mutant(sequence, position, new_base)
    wt = wild_type if (wild_type is not None and wild_type.sequence == sequence
                       and wild_type.algorithm == algorithm) else _fold(sequence, algorithm)
    mut = _fold(mutant_seq, algorithm)

    comparison = compare_pairs(wt.pairs, mut.pairs, len(sequence))
    wt_e, mut_e = _energy(wt), _energy(mut)
    ddg = round(mut_e - wt_e, 2) + 0.0
    if algorithm == "Thermodynamic":
        note = "Minimum estimated free energies from the simplified thermodynamic model."
    else:
        note = ("Nussinov does not use energies. These values score the two Nussinov "
                "structures with the thermodynamic energy model, for comparison only.")

    return MutationAnalysis(
        algorithm=algorithm,
        mutation=mutation,
        wild_type=wt,
        mutant=mut,
        comparison=comparison,
        wt_energy=wt_e,
        mut_energy=mut_e,
        ddg=ddg,
        energy_note=note,
        wt_elements=analyse_structure(sequence, wt.pairs).counts(),
        mut_elements=analyse_structure(mutant_seq, mut.pairs).counts(),
        wt_site=describe_site(wt, mutation.index),
        mut_site=describe_site(mut, mutation.index),
    )


# ------------------------------------------------------------------ #
# Interpretation (rule-based, only from computed values)
# ------------------------------------------------------------------ #
def interpret(a: MutationAnalysis) -> List[str]:
    m = a.mutation
    lines: List[str] = []
    wt_partner = pair_partner_map(a.wild_type.pairs, len(a.wild_type.sequence))[m.index]

    # 1. What happened at the mutated site
    if wt_partner >= 0:
        other = a.wild_type.sequence[wt_partner]
        if can_pair(m.replacement, other):
            lines.append(f"In the wild-type prediction, {m.original}{m.position} is {a.wt_site}. "
                         f"{m.replacement} can still pair with {other}{wt_partner + 1} "
                         f"({pair_type(m.replacement, other)} pair), so this pair is not "
                         "necessarily lost.")
        else:
            b1, b2 = (m.replacement, other) if m.index < wt_partner else (other, m.replacement)
            lines.append(f"In the wild-type prediction, {m.original}{m.position} is {a.wt_site}. "
                         f"{b1}-{b2} is not an allowed pair, so this base pair "
                         "cannot exist in the mutant.")
    else:
        lines.append(f"In the wild-type prediction, {m.original}{m.position} is {a.wt_site}. "
                     "Mutating an unpaired base does not directly remove a pair, but it can "
                     "allow new pairs or change loop energies.")

    # 2. Structural change
    c = a.comparison
    if c.identical:
        lines.append("The predicted mutant structure has exactly the same base pairs as the "
                     "wild type.")
    else:
        lines.append(f"The predicted structure changes: {len(c.retained)} pair(s) are retained, "
                     f"{len(c.lost)} lost and {len(c.new)} newly formed; the pairing status of "
                     f"{len(c.changed_positions)} nucleotide(s) differs.")
        if len(c.lost) + len(c.new) > 2:
            lines.append("Because more pairs change than the single mutated position can "
                         "explain, the fold has rearranged rather than only losing the pair at "
                         "the mutation site.")
    changed_elements = [k for k in a.wt_elements if a.wt_elements[k] != a.mut_elements.get(k)]
    if changed_elements:
        details = ", ".join(f"{k.lower()} {a.wt_elements[k]} → {a.mut_elements[k]}"
                            for k in changed_elements)
        lines.append(f"Structural elements changed: {details}.")

    # 3. Energy
    if abs(a.ddg) < ZERO_TOLERANCE:
        lines.append("ΔΔG = 0.00 kcal/mol: no change in the estimated folding free energy.")
    elif a.ddg > 0:
        lines.append(f"ΔΔG = {a.ddg:+.2f} kcal/mol: the mutant has a higher estimated folding "
                     "free energy, i.e. it is predicted to be less stable under this model.")
    else:
        lines.append(f"ΔΔG = {a.ddg:+.2f} kcal/mol: the mutant has a lower estimated folding "
                     "free energy, i.e. it is predicted to be more stable under this model.")
    if a.algorithm == "Nussinov":
        lines.append("Note: Nussinov maximises pair count, and ties are common, so some "
                     "differences may reflect an equally good alternative structure rather "
                     "than a real effect of the mutation.")

    # 4. Caveat
    lines.append("These are computational estimates from a simplified model. A change in "
                 "predicted structure or stability does not by itself show a loss or gain of "
                 "biological function.")
    return lines


# ------------------------------------------------------------------ #
# Mutation scan (optional extension)
# ------------------------------------------------------------------ #
def scan_mutations(sequence: str, algorithm: str = "Thermodynamic",
                   progress: Optional[Callable[[float], None]] = None) -> List[dict]:
    """
    Fold every possible single-nucleotide mutant (3 per position).
    Returns one row per mutant with ddG and pair changes.
    `progress` (optional) is called with a value from 0 to 1.
    """
    n = len(sequence)
    if n > SCAN_MAX_LENGTH:
        raise ValueError(f"The mutation scan is limited to {SCAN_MAX_LENGTH} nt "
                         f"({3 * n} folds would be needed for {n} nt).")
    wt = _fold(sequence, algorithm)
    wt_e = _energy(wt)
    wt_partner = pair_partner_map(wt.pairs, n)
    rows = []
    total, done = 3 * n, 0
    for pos in range(1, n + 1):
        for base in BASES:
            if base == sequence[pos - 1]:
                continue
            mutant_seq, mutation = create_mutant(sequence, pos, base)
            mut = _fold(mutant_seq, algorithm)
            comp = compare_pairs(wt.pairs, mut.pairs, n)
            mut_e = _energy(mut)
            rows.append({
                "Position": pos,
                "Original": mutation.original,
                "Replacement": base,
                "Mutation": mutation.label,
                "WT site paired": wt_partner[pos - 1] >= 0,
                "Pairs lost": len(comp.lost),
                "Pairs gained": len(comp.new),
                "Structure changed": not comp.identical,
                "Mutant ΔG (kcal/mol)": mut_e,
                "ΔΔG (kcal/mol)": round(mut_e - wt_e, 2) + 0.0,
            })
            done += 1
            if progress is not None:
                progress(done / total)
    return rows
