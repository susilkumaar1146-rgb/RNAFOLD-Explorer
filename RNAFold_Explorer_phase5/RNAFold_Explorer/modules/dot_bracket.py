"""
dot_bracket.py
==============
Conversion between a base-pair list and dot-bracket notation.

Dot-bracket notation
--------------------
    GGGGAAAACCCC
    ((((....))))
'(' marks the 5' partner of a pair, ')' the 3' partner, '.' an unpaired
base. Because only NESTED structures are allowed (no pseudoknots), one
bracket type is enough: a ')' always closes the most recent unmatched '('.
This is exactly how a stack works, which is how we parse it.

Functions
---------
pairs_to_dot_bracket(pairs, n)   -> str
dot_bracket_to_pairs(structure)  -> list of (i, j), 0-based
validate_dot_bracket(...)        -> list of error messages ([] = valid)
validate_pairs(...)              -> list of error messages for a pair list
"""

from typing import List, Optional, Tuple

from modules.pairing import MIN_HAIRPIN, can_pair

Pair = Tuple[int, int]


def validate_pairs(pairs: List[Pair], n: int, sequence: Optional[str] = None,
                   check_min_hairpin: bool = True) -> List[str]:
    """
    Check that a pair list describes a valid nested secondary structure:
      * indices are inside the sequence and i < j,
      * no position is used twice,
      * no two pairs cross (pseudoknot),
      * (optional) bases are complementary and hairpins have >= 3 nt.
    """
    errors: List[str] = []
    used = {}
    for i, j in pairs:
        if not (0 <= i < j < n):
            errors.append(f"Invalid pair ({i + 1}, {j + 1}) for length {n}.")
            continue
        for pos in (i, j):
            if pos in used:
                errors.append(f"Position {pos + 1} appears in more than one pair.")
            used[pos] = (i, j)
        if sequence is not None and not can_pair(sequence[i], sequence[j]):
            errors.append(f"Pair ({i + 1}, {j + 1}) {sequence[i]}-{sequence[j]} "
                          "is not A-U, G-C or G-U.")
        if check_min_hairpin and j - i - 1 < MIN_HAIRPIN:
            errors.append(f"Pair ({i + 1}, {j + 1}) encloses fewer than "
                          f"{MIN_HAIRPIN} nucleotides.")

    # Crossing check: pairs (i,j) and (k,l) cross if i < k < j < l.
    sorted_pairs = sorted(p for p in pairs if 0 <= p[0] < p[1] < n)
    for a in range(len(sorted_pairs)):
        i, j = sorted_pairs[a]
        for b in range(a + 1, len(sorted_pairs)):
            k, l = sorted_pairs[b]
            if k > j:
                break
            if i < k < j < l:
                errors.append(f"Pairs ({i + 1}, {j + 1}) and ({k + 1}, {l + 1}) "
                              "cross (pseudoknot); only nested structures are allowed.")
    return errors


def pairs_to_dot_bracket(pairs: List[Pair], n: int) -> str:
    """Build the dot-bracket string of length n from a (0-based) pair list."""
    errors = validate_pairs(pairs, n, check_min_hairpin=False)
    if errors:
        raise ValueError("Cannot build dot-bracket: " + "; ".join(errors))
    chars = ["."] * n
    for i, j in pairs:
        chars[i] = "("
        chars[j] = ")"
    return "".join(chars)


def dot_bracket_to_pairs(structure: str) -> List[Pair]:
    """
    Parse dot-bracket into a sorted list of 0-based pairs using a stack.
    Raises ValueError on unbalanced brackets or unknown characters.
    """
    stack: List[int] = []
    pairs: List[Pair] = []
    for pos, ch in enumerate(structure):
        if ch == "(":
            stack.append(pos)
        elif ch == ")":
            if not stack:
                raise ValueError(f"Unmatched ')' at position {pos + 1}.")
            pairs.append((stack.pop(), pos))
        elif ch != ".":
            raise ValueError(f"Invalid character '{ch}' at position {pos + 1}; "
                             "only '(', ')' and '.' are allowed.")
    if stack:
        raise ValueError(f"Unmatched '(' at position {stack[-1] + 1}.")
    return sorted(pairs)


def validate_dot_bracket(structure: str, sequence: Optional[str] = None,
                         check_min_hairpin: bool = True) -> List[str]:
    """Return a list of problems with a dot-bracket string ([] means valid)."""
    errors: List[str] = []
    if sequence is not None and len(structure) != len(sequence):
        errors.append(f"Structure length ({len(structure)}) differs from "
                      f"sequence length ({len(sequence)}).")
        return errors
    try:
        pairs = dot_bracket_to_pairs(structure)
    except ValueError as exc:
        return [str(exc)]
    errors.extend(validate_pairs(pairs, len(structure), sequence, check_min_hairpin))
    return errors


def pair_partner_map(pairs: List[Pair], n: int) -> List[int]:
    """partner[k] = index paired with k, or -1 if k is unpaired."""
    partner = [-1] * n
    for i, j in pairs:
        partner[i] = j
        partner[j] = i
    return partner
