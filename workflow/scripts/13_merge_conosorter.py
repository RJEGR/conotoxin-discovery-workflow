#!/usr/bin/env python3
# =============================================================================
# 26_stage7_merge_conosorter_simple.py  -  Stage 7 R4b: merge ConoSorter
#                                            Regex + pHMM tables into a single
#                                            non-redundant candidate set.
# =============================================================================
# Simplified replacement for legacy `2_merge_conosorter_results.py`. The
# legacy script merged 6 tables (3 assemblers × 2 modes) because ConoSorter
# ran per-assembler. In the new pipeline (agreed 2026-07-24), MMseqs2 dedup
# happens BEFORE ConoSorter, so there are only 2 tables per fold: one Regex
# and one pHMM run on the same dereplicated peptide set. That collapses the
# merge to a trivial join keyed on "Protein Sequence".
#
# INPUTS (via CLI):
#   --regex-tab   /path/to/Regex.tab   (from R4a)
#   --phmm-tab    /path/to/pHMM.tab    (from R4a)
#
# OUTPUTS (via CLI):
#   --out-tab     /path/to/merge.tab      : pipe-delimited, 1 row per unique
#                                            Protein Sequence, cols prefixed
#                                            with `regex__` and `phmm__`.
#   --out-fasta   /path/to/merge.fasta    : one FASTA per unique Protein Sequence
#                                            with synthetic ID `>seq_<idx>`.
#   --audit       /path/to/audit.tsv      : counts (n_regex, n_phmm, n_both,
#                                            n_regex_only, n_phmm_only,
#                                            n_unique_seqs).
#
# CANONICAL SCHEMA of ConoSorter output tables (both Regex and pHMM):
#   Conopeptide | Read Name | Sequence Frequency | Protein Sequence | # A.A |
#   # Cysteine(s) | % Hydrophobicity (Signal) | Superfamily (Signal) |
#   Superfamily (Pro-region) | Superfamily (Mature) | Score Superfamily |
#   Class (Signal) | Class (Pro-region) | Class (Mature) | Score Class
#
# The merged table keeps the shared cols (Protein Sequence, # A.A,
# # Cysteine(s), % Hydrophobicity, Sequence Frequency) unprefixed, and
# prefixes the mode-specific cols with `regex__` and `phmm__`.
# =============================================================================
# NOTE: no `from __future__ import annotations` — see gotcha #1 in
# snakemake_script_gotchas.md.
# =============================================================================

import argparse
import csv
import sys
from pathlib import Path


SEP = "|"

# Cols that are identical between Regex and pHMM tabs (same input seq) — kept
# unprefixed in the merged table.
SHARED_COLS = ["Protein Sequence", "# A.A", "# Cysteine(s)",
               "% Hydrophobicity (Signal)", "Sequence Frequency"]

# Cols that are mode-specific — prefixed `regex__` / `phmm__` in merged.
SPECIFIC_COLS = ["Conopeptide", "Read Name",
                 "Superfamily (Signal)", "Superfamily (Pro-region)",
                 "Superfamily (Mature)", "Score Superfamily",
                 "Class (Signal)", "Class (Pro-region)", "Class (Mature)",
                 "Score Class"]


def read_conosorter_tab(path: Path):
    """Return {protein_sequence: row_dict}. Handles empty files gracefully."""
    if not path.exists() or path.stat().st_size == 0:
        return {}, []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        first = fh.readline().rstrip("\n\r")
        if not first:
            return {}, []
        headers = [h.strip() for h in first.split(SEP)]
        rows_by_seq = {}
        for line_no, raw in enumerate(fh, start=2):
            line = raw.rstrip("\n\r")
            if not line:
                continue
            parts = line.split(SEP)
            # Pad or trim to header length
            if len(parts) < len(headers):
                parts = parts + [""] * (len(headers) - len(parts))
            elif len(parts) > len(headers):
                parts = parts[: len(headers) - 1] + [SEP.join(parts[len(headers) - 1:])]
            row = {h: v.strip() for h, v in zip(headers, parts)}
            seq = row.get("Protein Sequence", "").strip().upper()
            if not seq:
                continue
            # Dedup within a single tab by exact Protein Sequence
            if seq not in rows_by_seq:
                rows_by_seq[seq] = row
    return rows_by_seq, headers


def build_merged_header():
    cols = ["seq_id"] + SHARED_COLS[:]
    for m in ("regex", "phmm"):
        for c in SPECIFIC_COLS:
            cols.append(f"{m}__{c}")
    return cols


def build_merged_row(seq: str, idx: int, r_row: dict, p_row: dict) -> dict:
    """r_row / p_row may be None (only one mode had this seq)."""
    out = {"seq_id": f"seq_{idx:06d}"}
    # Shared cols: prefer regex row, fall back to phmm row
    donor = r_row or p_row
    for c in SHARED_COLS:
        out[c] = (donor or {}).get(c, "")
    # Mode-specific cols
    for c in SPECIFIC_COLS:
        out[f"regex__{c}"] = (r_row or {}).get(c, "")
        out[f"phmm__{c}"]  = (p_row or {}).get(c, "")
    # Ensure Protein Sequence is the canonical value (uppercased)
    out["Protein Sequence"] = seq
    return out


def write_merged_tab(rows: list, headers: list, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        fh.write(SEP.join(headers) + "\n")
        for row in rows:
            fh.write(SEP.join(str(row.get(h, "")).replace(SEP, "_")
                              for h in headers) + "\n")


def write_merged_fasta(rows: list, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        for row in rows:
            fh.write(f">{row['seq_id']}\n{row['Protein Sequence']}\n")


def write_audit(counts: dict, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("metric\tvalue\n")
        for k, v in counts.items():
            fh.write(f"{k}\t{v}\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--regex-tab", required=True, type=Path)
    ap.add_argument("--phmm-tab",  required=True, type=Path)
    ap.add_argument("--out-tab",   required=True, type=Path)
    ap.add_argument("--out-fasta", required=True, type=Path)
    ap.add_argument("--audit",     required=True, type=Path)
    args = ap.parse_args()

    r_rows, _ = read_conosorter_tab(args.regex_tab)
    p_rows, _ = read_conosorter_tab(args.phmm_tab)

    r_keys = set(r_rows.keys())
    p_keys = set(p_rows.keys())
    all_keys = sorted(r_keys | p_keys)

    n_regex      = len(r_keys)
    n_phmm       = len(p_keys)
    n_both       = len(r_keys & p_keys)
    n_regex_only = len(r_keys - p_keys)
    n_phmm_only  = len(p_keys - r_keys)
    n_unique     = len(all_keys)

    headers = build_merged_header()
    merged_rows = [
        build_merged_row(seq, i + 1, r_rows.get(seq), p_rows.get(seq))
        for i, seq in enumerate(all_keys)
    ]

    write_merged_tab(merged_rows, headers, args.out_tab)
    write_merged_fasta(merged_rows, args.out_fasta)
    write_audit({
        "n_regex_hits":      n_regex,
        "n_phmm_hits":       n_phmm,
        "n_regex_and_phmm":  n_both,
        "n_regex_only":      n_regex_only,
        "n_phmm_only":       n_phmm_only,
        "n_unique_seqs":     n_unique,
    }, args.audit)

    print(f"[merge_conosorter] n_unique_seqs={n_unique} "
          f"regex={n_regex} phmm={n_phmm} both={n_both}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
