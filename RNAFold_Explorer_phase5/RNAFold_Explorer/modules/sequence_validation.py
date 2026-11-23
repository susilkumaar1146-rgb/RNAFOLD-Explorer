"""
sequence_validation.py
======================
Cleans and validates user-supplied RNA sequences before any prediction.

Biology
-------
RNA contains A, U, G and C. DNA uses T instead of U, so a DNA-style
sequence can be analysed by converting every T to U (this models the
RNA transcript of that DNA sequence).

Rules implemented
-----------------
1. Upper-case all letters.
2. Remove whitespace (spaces, tabs, line breaks) - these are formatting,
   not sequence.
3. Convert T -> U (reported to the user, never silent).
4. Any other character (N, digits, '-', '*', ...) is an ERROR. We never
   silently drop characters, because that would shift every downstream
   position and make the base-pair numbers wrong.
5. The sequence must be at least MIN_LENGTH nt long: the smallest
   possible structure is one pair enclosing a 3-nt hairpin (1+3+1 = 5).

FASTA
-----
FASTA files contain one or more records:
    >header line
    SEQUENCE LINES...
parse_fasta() returns a list of (header, raw_sequence) tuples.

Output
------
validate_sequence() returns a ValidationResult object containing the
cleaned sequence, a validity flag, error messages and warnings.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from modules.pairing import MIN_HAIRPIN

VALID_BASES = set("AUGC")
MIN_LENGTH = MIN_HAIRPIN + 2  # 5 nt

# Explanations for commonly seen invalid characters.
_KNOWN_SYMBOLS = {
    "N": "ambiguous base (N)",
    "R": "IUPAC ambiguity code (A/G)",
    "Y": "IUPAC ambiguity code (C/U)",
    "-": "alignment gap",
    ".": "alignment gap or dot",
    "*": "stop/terminator symbol",
}


@dataclass
class ValidationResult:
    """Outcome of validating one sequence."""
    is_valid: bool
    sequence: str                      # cleaned sequence ('' if invalid)
    original_length: int               # characters before whitespace removal
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    t_converted: int = 0               # number of T -> U conversions

    @property
    def length(self) -> int:
        return len(self.sequence)


def clean_sequence(raw: str) -> Tuple[str, int]:
    """
    Upper-case, remove whitespace and convert T to U.
    Returns (cleaned_string, number_of_T_converted).
    Invalid characters are kept so that validate_sequence() can report them.
    """
    no_space = "".join(raw.split()).upper()
    t_count = no_space.count("T")
    return no_space.replace("T", "U"), t_count


def find_invalid_characters(sequence: str) -> List[Tuple[int, str]]:
    """Return (1-based position, character) for every non-AUGC character."""
    return [(pos + 1, ch) for pos, ch in enumerate(sequence) if ch not in VALID_BASES]


def validate_sequence(raw: Optional[str], max_length: Optional[int] = None) -> ValidationResult:
    """
    Validate a raw sequence string.

    Parameters
    ----------
    raw : str
        The text typed or uploaded by the user (without a FASTA header).
    max_length : int, optional
        Upper limit on length (algorithms set this; None = no limit).
    """
    if raw is None:
        raw = ""
    cleaned, t_count = clean_sequence(raw)
    errors: List[str] = []
    warnings: List[str] = []

    if cleaned == "":
        errors.append("The sequence is empty. Please enter an RNA sequence (A, U, G, C).")
        return ValidationResult(False, "", len(raw), errors, warnings, 0)

    if cleaned.startswith(">"):
        errors.append("The input looks like FASTA. Please use the FASTA upload "
                      "option or remove the '>' header line.")
        return ValidationResult(False, "", len(raw), errors, warnings, t_count)

    invalid = find_invalid_characters(cleaned)
    if invalid:
        shown = invalid[:10]
        details = ", ".join(
            f"'{ch}' at position {pos}" + (f" ({_KNOWN_SYMBOLS[ch]})" if ch in _KNOWN_SYMBOLS else "")
            for pos, ch in shown
        )
        more = f" ... and {len(invalid) - 10} more" if len(invalid) > 10 else ""
        errors.append(f"Invalid character(s) found: {details}{more}. "
                      "Only A, U, G and C are allowed (T is converted to U).")

    if len(cleaned) < MIN_LENGTH:
        errors.append(f"The sequence has {len(cleaned)} nt. At least {MIN_LENGTH} nt are "
                      f"needed to form one base pair enclosing a {MIN_HAIRPIN}-nt hairpin loop.")

    if max_length is not None and len(cleaned) > max_length:
        errors.append(f"The sequence has {len(cleaned)} nt, which exceeds the limit of "
                      f"{max_length} nt for this analysis.")

    if t_count:
        warnings.append(f"{t_count} T (thymine) base(s) were converted to U (uracil).")

    if errors:
        return ValidationResult(False, "", len(raw), errors, warnings, t_count)
    return ValidationResult(True, cleaned, len(raw), errors, warnings, t_count)


def gc_content(sequence: str) -> float:
    """Percentage of G + C in a cleaned sequence (0.0 for empty input)."""
    if not sequence:
        return 0.0
    gc = sequence.count("G") + sequence.count("C")
    return 100.0 * gc / len(sequence)


def base_composition(sequence: str) -> dict:
    """Count of each base, e.g. {'A': 4, 'U': 0, 'G': 4, 'C': 4}."""
    return {base: sequence.count(base) for base in "AUGC"}


def parse_fasta(text: str) -> List[Tuple[str, str]]:
    """
    Parse FASTA text into a list of (header, raw_sequence) records.

    * Lines starting with '>' start a new record.
    * Lines starting with ';' are comments and ignored.
    * Sequence lines are joined; validation happens later.
    * Text without any '>' header is treated as a single unnamed record.
    """
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    seq_lines: List[str] = []

    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_lines)))
            header = line[1:].strip() or f"sequence_{len(records) + 1}"
            seq_lines = []
        else:
            if header is None:
                header = "unnamed_sequence"
            seq_lines.append(line)

    if header is not None:
        records.append((header, "".join(seq_lines)))
    return records
