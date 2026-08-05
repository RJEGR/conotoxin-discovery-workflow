#!/usr/bin/env python3
# =============================================================================
# 32_stage7_classify_confidence.py  -  Stage 7 R9: confidence-tier
#                                       classification of annotated candidates.
# =============================================================================
# Thin I/O wrapper around the VENDORED legacy decision tree
# (`lib/legacy12_classify.py`, copied verbatim from legacy script 12 v7).
# The tree itself -- peptide-class-to-gene-superfamily translation, project
# direct equivalences, compatible-pair resolution, conditional-MEDIA relations
# and hard-conflict relations -- lives there untouched.
#
# TIERS (as described in the Stage 7 methods):
#   ALTA           ConoSorter-assigned AND ConoSorter agrees with >= 1 external
#                  source on the same superfamily.
#   MEDIA          Partial support, or consensus among external databases
#                  despite disagreement with ConoSorter.
#   BAJA           Only weak or conflicting support.
#   SIN_CLASIFICAR No meaningful agreement (includes rows whose labels failed
#                  normalization -- confirmed policy 2026-07-29).
#
# `is_positive` is taken from the `part` column: "positivo" == the
# ConoSorter-assigned set, "negativo" == the ConoSorter-unassigned set. The
# tree grades the two sets by different branches, exactly as the methods state
# ("ConoSorter-unassigned sequences were graded solely from similarity-based
# annotations").
#
# INPUTS  (via CLI): --normalized normalized.tsv          (from R8)
# OUTPUTS (via CLI): --out-xlsx  classification.xlsx      (4 tier sheets + all)
#                    --out-tsv   classification.tsv       (tidy, input to R10)
#                    --out-summary tier_summary.tsv
#                    --out-audit classify_audit.tsv
# =============================================================================
# NOTE: no `from __future__ import annotations` — see gotcha #1.
# =============================================================================

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

import openpyxl
from openpyxl.styles import Alignment, Font, PatternFill

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.legacy12_classify import classify_sequence   # noqa: E402

TIERS = ["ALTA", "MEDIA", "BAJA", "SIN_CLASIFICAR"]

SHEET_NAMES = {
    "ALTA": "Certidumbre_ALTA",
    "MEDIA": "Certidumbre_MEDIA",
    "BAJA": "Certidumbre_BAJA",
    "SIN_CLASIFICAR": "Sin_clasificar",
}

TIER_FILL = {
    "ALTA":           "C6EFCE",   # green
    "MEDIA":          "FFEB9C",   # amber
    "BAJA":           "FFC7CE",   # red
    "SIN_CLASIFICAR": "D9D9D9",   # grey
}

OUT_COLS = [
    "seq_id", "part",
    "Certidumbre", "Regla de clasificacion", "Consenso",
    "Superfamily (Signal)", "conoserver", "ncbi", "uniprot",
    "raw__Superfamily (Signal)", "raw__conoserver", "raw__ncbi", "raw__uniprot",
    "Protein Sequence", "# A.A", "# Cysteine(s)",
]


def classify_rows(rows: list, audit: Counter):
    out = []
    for r in rows:
        is_positive = r.get("part", "") == "positivo"
        certidumbre, regla, consenso = classify_sequence(
            is_positive=is_positive,
            conosorter=r.get("Superfamily (Signal)", ""),
            conoserver=r.get("conoserver", ""),
            ncbi=r.get("ncbi", ""),
            uniprot=r.get("uniprot", ""),
        )
        audit[f"tier_{certidumbre}"] += 1
        audit[f"tier_{certidumbre}__{r.get('part', '?')}"] += 1
        # Rule families are useful for QC; keep only the stem before the
        # nomenclature-mapping suffix so the histogram stays readable.
        audit[f"rule_{regla.split('__mapeo_nomenclatural_')[0]}"] += 1
        d = dict(r)
        d["Certidumbre"] = certidumbre
        d["Regla de clasificacion"] = regla
        d["Consenso"] = consenso
        out.append(d)
    return out


def write_xlsx(rows: list, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    header_font = Font(bold=True, color="FFFFFF")
    header_fill = PatternFill("solid", fgColor="4472C4")

    def add_sheet(title: str, subset: list, fill_hex: str = None):
        ws = wb.create_sheet(title[:31])
        ws.append(OUT_COLS)
        for c in range(1, len(OUT_COLS) + 1):
            cell = ws.cell(row=1, column=c)
            cell.font = header_font
            cell.fill = header_fill
            cell.alignment = Alignment(vertical="center", wrap_text=True)
        for r in subset:
            ws.append([r.get(c, "") for c in OUT_COLS])
        if fill_hex and subset:
            f = PatternFill("solid", fgColor=fill_hex)
            for row in range(2, len(subset) + 2):
                ws.cell(row=row, column=3).fill = f
        ws.freeze_panes = "A2"
        for i, col in enumerate(OUT_COLS, start=1):
            ws.column_dimensions[
                openpyxl.utils.get_column_letter(i)].width = min(max(len(col) + 2, 12), 60)

    for tier in TIERS:
        add_sheet(SHEET_NAMES[tier],
                  [r for r in rows if r["Certidumbre"] == tier],
                  TIER_FILL[tier])
    add_sheet("Todas_las_secuencias", rows)
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
    ap.add_argument("--normalized",  required=True, type=Path)
    ap.add_argument("--out-xlsx",    required=True, type=Path)
    ap.add_argument("--out-tsv",     required=True, type=Path)
    ap.add_argument("--out-summary", required=True, type=Path)
    ap.add_argument("--out-audit",   required=True, type=Path)
    args = ap.parse_args()

    with args.normalized.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    audit = Counter()
    classified = classify_rows(rows, audit)

    write_xlsx(classified, args.out_xlsx)
    write_tsv(classified, OUT_COLS, args.out_tsv)

    # ---- Tier summary (tidy; feeds the Stage 7 manifest) -----------------
    summary = []
    for tier in TIERS:
        for part in ("positivo", "negativo"):
            n = sum(1 for r in classified
                    if r["Certidumbre"] == tier and r["part"] == part)
            summary.append({"tier": tier, "part": part, "n": n})
        summary.append({"tier": tier, "part": "all",
                        "n": sum(1 for r in classified
                                 if r["Certidumbre"] == tier)})
    write_tsv(summary, ["tier", "part", "n"], args.out_summary)

    args.out_audit.parent.mkdir(parents=True, exist_ok=True)
    with args.out_audit.open("w", encoding="utf-8") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"rows_total\t{len(classified)}\n")
        for k in sorted(audit):
            fh.write(f"{k}\t{audit[k]}\n")

    tiers = {t: audit.get(f"tier_{t}", 0) for t in TIERS}
    print(f"[classify_confidence] n={len(classified)} tiers={tiers}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
