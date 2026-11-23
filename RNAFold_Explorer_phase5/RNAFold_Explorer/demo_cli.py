"""
demo_cli.py
===========
Command-line demonstration of the Phase 2 algorithms (no Streamlit needed).

Usage (from the RNAFold_Explorer folder):
    python demo_cli.py                            # runs the example FASTA file
    python demo_cli.py GGGGAAAACCCC               # folds your own sequence
    python demo_cli.py --mutate GGGGAAAACCCC C12A # mutation analysis
"""

import os
import sys

from modules import nussinov, thermodynamic
from modules.energy_eval import evaluate_structure
from modules.mutation_analysis import analyse_mutation, interpret, parse_mutation
from modules.sequence_validation import gc_content, parse_fasta, validate_sequence
from modules.structure_elements import analyse_structure

HERE = os.path.dirname(os.path.abspath(__file__))


def report(name: str, raw: str) -> None:
    print("=" * 78)
    print(name)
    check = validate_sequence(raw, max_length=thermodynamic.MAX_LENGTH)
    for w in check.warnings:
        print("  warning:", w)
    if not check.is_valid:
        for e in check.errors:
            print("  ERROR:", e)
        return
    seq = check.sequence
    print(f"Length {len(seq)} nt, GC content {gc_content(seq):.1f}%")
    print("Sequence      ", seq)

    nus = nussinov.fold(seq)
    print("Nussinov      ", nus.dot_bracket, f"  {nus.pair_count} pairs")
    for note in nus.notes[1:]:
        print("               ", note)

    th = thermodynamic.fold(seq)
    print("Thermodynamic ", th.dot_bracket, f"  {th.energy:.2f} kcal/mol (estimated)")

    print("Loop-by-loop energy of the thermodynamic structure:")
    for loop in evaluate_structure(seq, th.pairs).loops:
        if loop.energy != 0 or loop.loop_type != "Exterior":
            print(f"   {loop.energy:+6.2f}  {loop.describe(seq)}")

    counts = analyse_structure(seq, th.pairs).counts()
    print("Structural elements:", ", ".join(f"{k}: {v}" for k, v in counts.items()))

    nus_energy = evaluate_structure(seq, nus.pairs).total
    print(f"Energy of the Nussinov structure under the same model: {nus_energy:.2f} kcal/mol")


def mutation_report(raw: str, mutation_text: str) -> None:
    check = validate_sequence(raw, max_length=thermodynamic.MAX_LENGTH)
    if not check.is_valid:
        print("\n".join("ERROR: " + e for e in check.errors))
        return
    seq = check.sequence
    try:
        m = parse_mutation(mutation_text, seq)
    except ValueError as exc:
        print("ERROR:", exc)
        return
    for algorithm in ("Nussinov", "Thermodynamic"):
        a = analyse_mutation(seq, m.position, m.replacement, algorithm)
        c = a.comparison
        print("=" * 78)
        print(f"{algorithm}: mutation {m.label}")
        print("Wild type ", seq)
        print("          ", a.wild_type.dot_bracket, f"  {a.wild_type.pair_count} pairs")
        print("Mutant    ", a.mutant.sequence)
        print("          ", a.mutant.dot_bracket, f"  {a.mutant.pair_count} pairs")
        marker = "".join("^" if k == m.index else ("*" if k in c.changed_positions else " ")
                         for k in range(len(seq)))
        print("          ", marker, "  (^ mutation, * pairing changed)")
        print(f"Pairs retained {len(c.retained)}, lost {len(c.lost)}, new {len(c.new)}")
        print(f"Estimated dG: WT {a.wt_energy:.2f}, mutant {a.mut_energy:.2f}, "
              f"ddG {a.ddg:+.2f} kcal/mol")
        for line in interpret(a):
            print("  -", line)


def main() -> None:
    if len(sys.argv) == 4 and sys.argv[1] == "--mutate":
        mutation_report(sys.argv[2], sys.argv[3])
        return
    if len(sys.argv) > 1:
        report("User sequence", "".join(sys.argv[1:]))
        return
    with open(os.path.join(HERE, "data", "sample_sequences.fasta"), encoding="utf-8") as fh:
        for header, raw in parse_fasta(fh.read()):
            report(header, raw)


if __name__ == "__main__":
    main()
