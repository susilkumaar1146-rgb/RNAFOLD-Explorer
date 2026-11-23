"""
results_export.py
=================
Builds the downloadable files (as text or bytes) and the aligned text blocks
shown on screen. Pure Python: the Streamlit pages only pass these strings to
st.download_button / st.code.

Formats
-------
* FASTA            >name / sequence wrapped at 60 nt
* Dot-bracket      Vienna-style: >name, sequence, structure ( energy)
* CSV              base-pair tables, mutation comparison, mutation scan
* Text report      everything in one readable summary (with model caveats)
* ZIP              all of the above in one archive
"""

import csv
import io
import zipfile
from datetime import datetime
from typing import Dict, List, Optional

from modules.fold_result import FoldResult
from modules.sequence_validation import base_composition, gc_content
from modules.structure_elements import analyse_structure

DISCLAIMER = ("All structures and energies are computational estimates (Nussinov "
              "maximum-pair algorithm and a simplified nearest-neighbour model based on a subset "
              "of Turner 2004 parameters). They are not experimentally determined structures, "
              "and pseudoknots are not predicted.")


def safe_name(name: str) -> str:
    """Turn a sequence name into a short, file-system-safe file name."""
    keep = "".join(c if c.isalnum() or c in "-_." else "_" for c in name.strip())
    while "__" in keep:
        keep = keep.replace("__", "_")
    return keep[:40].strip("_.") or "sequence"


def fasta(name: str, sequence: str, width: int = 60) -> str:
    lines = [f">{name}"]
    lines += [sequence[k:k + width] for k in range(0, len(sequence), width)]
    return "\n".join(lines) + "\n"


def dot_bracket_text(name: str, result: FoldResult) -> str:
    """Vienna-style output (one line per item, energy in brackets when available)."""
    structure = result.dot_bracket
    if result.energy is not None:
        structure += f" ({result.energy:6.2f})"
    else:
        structure += f" [{result.pair_count} pairs]"
    return f">{name} | {result.algorithm}\n{result.sequence}\n{structure}\n"


def _csv(rows: List[dict]) -> str:
    if not rows:
        return ""
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0].keys()), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def pairs_csv(result: FoldResult) -> str:
    rows = result.pair_table()
    if not rows:
        return "Position i,Base i,Position j,Base j,Pair,Type\n"
    return _csv(rows)


def mutation_csv(analysis) -> str:
    """Summary values followed by the retained / lost / new pair lists."""
    rows = [{"Section": "Summary", "Item": r["Quantity"], "Value": r["Value"], "Unit": r["Unit"]}
            for r in analysis.summary_rows()]
    seq_wt, seq_mut = analysis.wild_type.sequence, analysis.mutant.sequence
    for label, pairs, seq in (("Retained pair", analysis.comparison.retained, seq_wt),
                              ("Lost pair", analysis.comparison.lost, seq_wt),
                              ("New pair", analysis.comparison.new, seq_mut)):
        for i, j in pairs:
            rows.append({"Section": label, "Item": f"{i + 1}-{j + 1}",
                         "Value": f"{seq[i]}-{seq[j]}", "Unit": ""})
    rows.append({"Section": "Structure", "Item": "Wild type", "Value": analysis.wild_type.dot_bracket, "Unit": ""})
    rows.append({"Section": "Structure", "Item": "Mutant", "Value": analysis.mutant.dot_bracket, "Unit": ""})
    return _csv(rows)


def scan_csv(rows: List[dict]) -> str:
    return _csv(rows)


def aligned_block(tracks: Dict[str, str], width: int = 60, ruler: bool = True) -> str:
    """
    Show several equal-length strings (sequence, structures, markers) one under
    another, wrapped every `width` columns with position numbers.
    Example:
                  1         11
        Sequence  GGGGAAAACCCC
        Structure ((((....))))
    """
    if not tracks:
        return ""
    length = len(next(iter(tracks.values())))
    label_width = max(len(k) for k in tracks) + 2
    out: List[str] = []
    for start in range(0, length, width):
        end = min(start + width, length)
        if ruler:
            marks = ""
            for pos in range(start + 1, end + 1):
                if pos == start + 1 or pos % 10 == 1:
                    marks = marks.ljust(pos - start - 1) + str(pos)
            out.append(" " * label_width + marks[: end - start + 6])
        for name, text in tracks.items():
            out.append(name.ljust(label_width) + text[start:end])
        out.append("")
    return "\n".join(out).rstrip() + "\n"


def summary_report(name: str, sequence: str, results: Dict[str, FoldResult],
                   mutation=None, comparison=None, scan_rows: Optional[List[dict]] = None) -> str:
    lines = [
        "RNAFold Explorer - prediction summary",
        "=" * 60,
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"Sequence name: {name}",
        f"Length: {len(sequence)} nt   GC content: {gc_content(sequence):.1f}%   "
        + "   ".join(f"{b}: {c}" for b, c in base_composition(sequence).items()),
        "",
        aligned_block({"Sequence": sequence}).rstrip(),
        "",
    ]
    for algorithm, result in results.items():
        lines.append(f"--- {algorithm} prediction ---")
        lines.append(aligned_block({"Sequence": sequence, "Structure": result.dot_bracket}).rstrip())
        lines.append(f"Base pairs: {result.pair_count}")
        if result.energy is not None:
            lines.append(f"Estimated free energy: {result.energy:.2f} kcal/mol")
        counts = analyse_structure(sequence, result.pairs).counts()
        lines.append("Structural elements: " + ", ".join(f"{k} {v}" for k, v in counts.items()))
        lines.append("Pairs (1-based): " + (", ".join(f"{i}-{j}" for i, j in result.pairs_1based()) or "none"))
        lines.extend(f"Note: {note}" for note in result.notes)
        lines.append("")
    if comparison is not None:
        lines.append("--- Algorithm comparison ---")
        lines.extend(comparison.explanation())
        lines.append("")
    if mutation is not None:
        a = mutation
        lines.append(f"--- Mutation analysis: {a.mutation.label} ({a.algorithm}) ---")
        lines.append(aligned_block({"Wild type": a.wild_type.sequence, "WT struct": a.wild_type.dot_bracket,
                                    "Mutant": a.mutant.sequence, "Mut struct": a.mutant.dot_bracket}).rstrip())
        for row in a.summary_rows():
            lines.append(f"{row['Quantity']}: {row['Value']} {row['Unit']}".rstrip())
        lines.append("Lost pairs: " + (", ".join(f"{i + 1}-{j + 1}" for i, j in a.comparison.lost) or "none"))
        lines.append("New pairs: " + (", ".join(f"{i + 1}-{j + 1}" for i, j in a.comparison.new) or "none"))
        lines.append(f"Energy note: {a.energy_note}")
        from modules.mutation_analysis import interpret
        lines.extend(f"* {line}" for line in interpret(a))
        lines.append("")
    if scan_rows:
        top = sorted(scan_rows, key=lambda r: r["ΔΔG (kcal/mol)"], reverse=True)[:10]
        lines.append("--- Mutation scan: 10 mutations with the largest ΔΔG ---")
        for r in top:
            lines.append(f"{r['Mutation']:>8}  ΔΔG {r['ΔΔG (kcal/mol)']:+.2f} kcal/mol   "
                         f"pairs lost {r['Pairs lost']}, gained {r['Pairs gained']}")
        lines.append("")
    lines.append("Disclaimer: " + DISCLAIMER)
    return "\n".join(lines) + "\n"


def bundle_zip(files: Dict[str, str]) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for filename, content in files.items():
            archive.writestr(filename, content)
    return buffer.getvalue()
