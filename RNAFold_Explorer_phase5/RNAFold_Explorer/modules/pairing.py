"""
pairing.py
==========
Base-pairing rules shared by every algorithm in RNAFold Explorer.

Biology
-------
RNA bases pair through hydrogen bonds:
    * A-U  (Watson-Crick, 2 hydrogen bonds)
    * G-C  (Watson-Crick, 3 hydrogen bonds)
    * G-U  (wobble pair, 2 hydrogen bonds, geometrically shifted)
A hairpin loop needs at least 3 unpaired nucleotides, otherwise the
backbone cannot bend back far enough for the two bases to meet.

Keeping these rules in ONE module guarantees that Nussinov, the
thermodynamic model, the energy evaluator and the validators all use
exactly the same definition of a "valid pair".
"""

# All six allowed ordered pairs (5' base, 3' base).
ALLOWED_PAIRS = {
    ("A", "U"), ("U", "A"),
    ("G", "C"), ("C", "G"),
    ("G", "U"), ("U", "G"),
}

# Pairs that receive the terminal AU/GU helix-end penalty in the energy model.
AU_GU_PAIRS = {("A", "U"), ("U", "A"), ("G", "U"), ("U", "G")}

# Minimum number of unpaired nucleotides enclosed by a hairpin.
MIN_HAIRPIN = 3


def can_pair(base1: str, base2: str) -> bool:
    """Return True if the two bases can form an A-U, G-C or G-U pair."""
    return (base1, base2) in ALLOWED_PAIRS


def can_pair_positions(sequence: str, i: int, j: int) -> bool:
    """
    Return True if positions i and j (0-based, i < j) may pair.

    Two conditions must hold:
      1. the bases are complementary (A-U, G-C or G-U), and
      2. at least MIN_HAIRPIN nucleotides lie between them.
    """
    if i < 0 or j >= len(sequence) or i >= j:
        return False
    if j - i - 1 < MIN_HAIRPIN:
        return False
    return can_pair(sequence[i], sequence[j])


def pair_type(base1: str, base2: str) -> str:
    """Classify a pair as 'Watson-Crick', 'Wobble' or 'Invalid'."""
    if (base1, base2) in {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}:
        return "Watson-Crick"
    if (base1, base2) in {("G", "U"), ("U", "G")}:
        return "Wobble"
    return "Invalid"


def is_au_gu(base1: str, base2: str) -> bool:
    """True for A-U, U-A, G-U and U-G pairs (weaker helix ends)."""
    return (base1, base2) in AU_GU_PAIRS
