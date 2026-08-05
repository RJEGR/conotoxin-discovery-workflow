#!/usr/bin/env python3
# =============================================================================
# 21_check_invariants.py  -  Stage 02 R10b: final QC gate
# =============================================================================
# Verifies the pipeline's output contract and writes a PASS/FAIL verdict.
#
# ALWAYS EXITS 0. A companion `shell:` gate rule reads the verdict and halts
# the DAG, so the diagnostic report survives on disk instead of being deleted
# by Snakemake as a failed job's output.
#
# CHECKS
#   1. Manifest exists, is non-empty, and carries the expected columns.
#   2. Row count >= require_manifest_min_rows.
#   3. Per group: nr.pep is non-empty (an empty dereplicated peptide set means
#      nothing downstream can be right).
#   4. Per group: SignalP positives — configurable. OFF by default for real
#      data: a library from the wrong tissue, or a failed prep, can honestly
#      yield zero secreted candidates, and that should be a reported note
#      rather than a failed run.
#   5. Per group: tier counts sum to the normalized row count. This catches
#      silent row loss between R8 and R9, which no other check would notice.
#   6. ADVISORY: unresolved-label burden versus max_conflicting_labels_pct.
#
# The distinction between 3 and 4 is the point of this gate: 3 is an
# arithmetic impossibility, 4 is a biological outcome. Only the first should
# ever stop a run.
# =============================================================================
# NOTE: no `from __future__ import annotations` — see project gotchas.
# =============================================================================

import argparse
import csv
import sys
from pathlib import Path


REQUIRED_COLS = {
    "group", "n_nr_orfs", "n_candidates", "n_signalp_positive",
    "n_normalized_rows", "n_alta", "n_media", "n_baja", "n_sin_clasificar",
}


def to_int(v, default=0):
    try:
        return int(float(str(v).strip()))
    except (TypeError, ValueError):
        return default


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest",    required=True, type=Path)
    ap.add_argument("--out-report",  required=True, type=Path)
    ap.add_argument("--out-verdict", required=True, type=Path)
    ap.add_argument("--min-rows", type=int, default=1)
    ap.add_argument("--require-nonempty-nr-pep",  default="true")
    ap.add_argument("--require-nonempty-signalp", default="false")
    ap.add_argument("--max-conflicting-pct", type=float, default=25.0)
    args = ap.parse_args()

    req_nr = str(args.require_nonempty_nr_pep).lower() == "true"
    req_sp = str(args.require_nonempty_signalp).lower() == "true"

    findings, fails = [], 0

    def add(check, status, detail):
        nonlocal fails
        findings.append((check, status, detail))
        if status == "FAIL":
            fails += 1

    # ---- 1. manifest + schema ---------------------------------------------
    rows = []
    if not args.manifest.exists() or not args.manifest.stat().st_size:
        add("manifest_exists", "FAIL", f"missing or empty: {args.manifest}")
    else:
        add("manifest_exists", "PASS", f"{args.manifest.stat().st_size} bytes")
        with args.manifest.open(encoding="utf-8", newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            cols = set(rd.fieldnames or [])
            rows = list(rd)
        missing = REQUIRED_COLS - cols
        if missing:
            add("manifest_schema", "FAIL", f"missing columns: {sorted(missing)}")
        else:
            add("manifest_schema", "PASS", f"{len(cols)} columns")

    # ---- 2. row count ------------------------------------------------------
    if len(rows) >= args.min_rows:
        add("manifest_min_rows", "PASS", f"{len(rows)} >= {args.min_rows}")
    else:
        add("manifest_min_rows", "FAIL", f"{len(rows)} < {args.min_rows}")

    # ---- 3..6 per group ----------------------------------------------------
    for r in rows:
        g = r.get("group", "?")

        if req_nr:
            n = to_int(r.get("n_nr_orfs"))
            add(f"nr_pep[{g}]", "PASS" if n > 0 else "FAIL", f"n_nr_orfs={n}")

        n_sp = to_int(r.get("n_signalp_positive"))
        if req_sp:
            add(f"signalp[{g}]", "PASS" if n_sp > 0 else "FAIL",
                f"n_signalp_positive={n_sp}")
        elif n_sp == 0:
            add(f"signalp[{g}]", "WARN",
                "no SignalP-positive candidates — plausible for a non-venom "
                "tissue or a failed library, but verify before interpreting")
        else:
            add(f"signalp[{g}]", "PASS", f"n_signalp_positive={n_sp}")

        norm = to_int(r.get("n_normalized_rows"))
        tiers = sum(to_int(r.get(k)) for k in
                    ("n_alta", "n_media", "n_baja", "n_sin_clasificar"))
        if norm == tiers:
            add(f"tier_arithmetic[{g}]", "PASS",
                f"{tiers} tier rows == {norm} normalized rows")
        else:
            add(f"tier_arithmetic[{g}]", "FAIL",
                f"tier sum {tiers} != normalized rows {norm}")

        unres = to_int(r.get("n_unresolved_labels"))
        pct = (unres / norm * 100) if norm else 0.0
        add(f"unresolved_labels[{g}]",
            "ADVISORY" if pct <= args.max_conflicting_pct else "WARN",
            f"{unres}/{norm} = {pct:.1f}% (threshold {args.max_conflicting_pct}%)")

        if r.get("notes"):
            add(f"notes[{g}]", "ADVISORY", r["notes"])

    verdict = "PASS" if fails == 0 else "FAIL"

    args.out_report.parent.mkdir(parents=True, exist_ok=True)
    with args.out_report.open("w", encoding="utf-8") as fh:
        fh.write("# Stage 02 annotation invariants\n")
        fh.write(f"# manifest = {args.manifest}\n")
        fh.write(f"# groups   = {len(rows)}\n\n")
        w = max((len(c) for c, _s, _d in findings), default=20) + 2
        for check, status, detail in findings:
            fh.write(f"[{status:8}] {check:<{w}} :: {detail}\n")
        fh.write(f"\nTotal FAIL: {fails}\nVERDICT: {verdict}\n")

    args.out_verdict.parent.mkdir(parents=True, exist_ok=True)
    args.out_verdict.write_text(
        f"{verdict}\ntotal_fails={fails}\ngroups={len(rows)}\n")

    print(f"[check_invariants] verdict={verdict} fails={fails} -> "
          f"{args.out_report}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
