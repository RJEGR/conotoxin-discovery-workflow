#!/usr/bin/env python3
# =============================================================================
# 31_stage7_normalize_final.py  -  Stage 7 R8: normalize superfamily/class
#                                   labels to a common vocabulary.
# =============================================================================
# Thin I/O wrapper around the VENDORED legacy normalization layer
# (`lib/legacy11_normalize.py`, copied verbatim from legacy script 11 v6).
# All curated exact mappings and rule-based transformations live there; this
# file only reads the annotated workbooks, applies the four normalizers, and
# writes the normalized outputs.
#
# WHY NORMALIZATION IS REQUIRED (per the Stage 7 methods): "Differences in
# notation, spacing, capitalization, and label wording produced apparent
# discrepancies among sources." Concretely, for the same candidate the four
# evidence sources may report:
#     ConoSorter Superfamily (Signal) : "A"
#     ConoServer                      : "A superfamily"
#     NCBI                            : "" (feature table rarely states it)
#     UniProt                         : "Alpha-conotoxin CIA (Cl)"
# All four must collapse to the token "A" before R9 can measure agreement.
#
# COLUMN DERIVATION. Legacy operated on a single `Superfamily (Signal)` column.
# The new pipeline runs ConoSorter in both modes, so that column is derived
# here: RegEx value first, pHMM as fallback (RegEx is the more sensitive mode
# and was the legacy default; the choice is recorded in the audit).
#
# UNRESOLVED LABELS. Values that survive normalization without landing in
# VALID_FINAL are left untouched and reported in `unresolved_labels.tsv`.
# They flow into R9 where the decision tree routes them to the
# `SIN_CLASIFICAR` tier. Confirmed policy 2026-07-29: do NOT invent mappings.
#
# INPUTS (via CLI):
#   --pos-xlsx / --neg-xlsx   annotate/{positivo,negativo}_annotated.xlsx
# OUTPUTS (via CLI):
#   --out-pos-xlsx / --out-neg-xlsx   normalized workbooks
#   --out-normalized                  normalized.tsv (tidy, input to R9)
#   --out-unresolved                  unresolved_labels.tsv
#   --out-audit                       normalize_audit.tsv
# =============================================================================
# NOTE: no `from __future__ import annotations` — see gotcha #1.
# =============================================================================

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

import openpyxl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.legacy11_normalize import (            # noqa: E402
    VALID_FINAL,
    clean_conoserver_value,
    clean_signal_value,
    normalize_ncbi_value,
    normalize_uniprot_v3,
)

KEY_COL = "seq_id"
SIGNAL_OUT_COL = "Superfamily (Signal)"
NORM_COLS = [SIGNAL_OUT_COL, "conoserver", "ncbi", "uniprot"]

# Source columns emitted by R4b / R7 in the workbooks.
SIGNAL_SRC_REGEX = "regex__Superfamily (Signal)"
SIGNAL_SRC_PHMM = "phmm__Superfamily (Signal)"

EMPTYISH = {"", "-", "na", "n/a", "none", "null", "nan"}


def _is_empty(v) -> bool:
    return str(v or "").strip().lower() in EMPTYISH


def read_principal(xlsx: Path):
    """Return (headers, rows) from the 'principal' sheet (or the first one)."""
    wb = openpyxl.load_workbook(xlsx, data_only=True)
    ws = wb["principal"] if "principal" in wb.sheetnames else wb[wb.sheetnames[0]]
    rows_iter = ws.iter_rows(values_only=True)
    headers = [str(h).strip() if h is not None else "" for h in next(rows_iter, [])]
    rows = []
    for r in rows_iter:
        if r is None:
            continue
        d = {h: ("" if v is None else str(v).strip())
             for h, v in zip(headers, r)}
        if d.get(KEY_COL):
            rows.append(d)
    wb.close()
    return headers, rows


def derive_signal(row: dict, audit: Counter) -> str:
    """RegEx value first, pHMM as fallback. Records which mode supplied it."""
    rx = row.get(SIGNAL_SRC_REGEX, "")
    if not _is_empty(rx):
        audit["signal_from_regex"] += 1
        return rx
    ph = row.get(SIGNAL_SRC_PHMM, "")
    if not _is_empty(ph):
        audit["signal_from_phmm"] += 1
        return ph
    audit["signal_absent"] += 1
    return ""


def normalize_row(row: dict, part: str, audit: Counter, unresolved: list):
    """Apply the four vendored normalizers. Returns the normalized dict."""
    raw_signal = derive_signal(row, audit)
    raw_cs = row.get("conoserver", "")
    raw_nc = row.get("ncbi", "")
    raw_up = row.get("uniprot", "")

    n_signal = clean_signal_value(raw_signal) if not _is_empty(raw_signal) else ""
    n_cs = clean_conoserver_value(raw_cs) if not _is_empty(raw_cs) else ""
    n_nc, st_nc = normalize_ncbi_value(raw_nc) if not _is_empty(raw_nc) else ("", "empty")
    n_up, st_up = normalize_uniprot_v3(raw_up) if not _is_empty(raw_up) else ("", "empty")

    out = {
        KEY_COL: row.get(KEY_COL, ""),
        "part": part,
        SIGNAL_OUT_COL: n_signal,
        "conoserver": n_cs,
        "ncbi": n_nc,
        "uniprot": n_up,
        f"raw__{SIGNAL_OUT_COL}": raw_signal,
        "raw__conoserver": raw_cs,
        "raw__ncbi": raw_nc,
        "raw__uniprot": raw_up,
        "status__ncbi": st_nc,
        "status__uniprot": st_up,
        "Protein Sequence": row.get("Protein Sequence", ""),
        "# A.A": row.get("# A.A", ""),
        "# Cysteine(s)": row.get("# Cysteine(s)", ""),
    }

    for col, raw in ((SIGNAL_OUT_COL, raw_signal), ("conoserver", raw_cs),
                     ("ncbi", raw_nc), ("uniprot", raw_up)):
        val = out[col]
        if not val:
            audit[f"{col}__empty"] += 1
        elif val in VALID_FINAL:
            audit[f"{col}__resolved"] += 1
        else:
            audit[f"{col}__unresolved"] += 1
            unresolved.append({"seq_id": out[KEY_COL], "part": part,
                               "source": col, "raw": raw, "normalized": val})
    return out


OUT_COLS = [KEY_COL, "part"] + NORM_COLS + [
    f"raw__{SIGNAL_OUT_COL}", "raw__conoserver", "raw__ncbi", "raw__uniprot",
    "status__ncbi", "status__uniprot",
    "Protein Sequence", "# A.A", "# Cysteine(s)",
]


def write_xlsx(rows: list, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "principal"
    ws.append(OUT_COLS)
    for r in rows:
        ws.append([r.get(c, "") for c in OUT_COLS])
    wb.save(path)


def write_tsv(rows: list, cols: list, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pos-xlsx",       required=True, type=Path)
    ap.add_argument("--neg-xlsx",       required=True, type=Path)
    ap.add_argument("--out-pos-xlsx",   required=True, type=Path)
    ap.add_argument("--out-neg-xlsx",   required=True, type=Path)
    ap.add_argument("--out-normalized", required=True, type=Path)
    ap.add_argument("--out-unresolved", required=True, type=Path)
    ap.add_argument("--out-audit",      required=True, type=Path)
    args = ap.parse_args()

    audit = Counter()
    unresolved = []
    all_rows = []

    for part, src, dst in (("positivo", args.pos_xlsx, args.out_pos_xlsx),
                           ("negativo", args.neg_xlsx, args.out_neg_xlsx)):
        _headers, rows = read_principal(src)
        norm = [normalize_row(r, part, audit, unresolved) for r in rows]
        audit[f"rows_{part}"] = len(norm)
        write_xlsx(norm, dst)
        all_rows.extend(norm)

    write_tsv(all_rows, OUT_COLS, args.out_normalized)
    write_tsv(unresolved,
              ["seq_id", "part", "source", "raw", "normalized"],
              args.out_unresolved)

    args.out_audit.parent.mkdir(parents=True, exist_ok=True)
    with args.out_audit.open("w", encoding="utf-8") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"rows_total\t{len(all_rows)}\n")
        for k in sorted(audit):
            fh.write(f"{k}\t{audit[k]}\n")
        fh.write(f"unresolved_rows\t{len(unresolved)}\n")
        fh.write(f"unresolved_distinct\t"
                 f"{len({(u['source'], u['raw']) for u in unresolved})}\n")

    print(f"[normalize_final] rows={len(all_rows)} unresolved={len(unresolved)}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
