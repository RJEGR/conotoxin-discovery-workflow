#!/usr/bin/env python3
# =============================================================================
# 20_aggregate_manifest.py  -  Stage 02 R10a: annotation_manifest.tsv
# =============================================================================
# One row per GROUP: the whole annotation audit trail for an assembly unit,
# flattened into the stage's output contract.
#
# The row key is the group because that is what was annotated: Stage 01
# assembled the group's pooled reads with each assembler, Stage 02 dereplicated
# those assemblies into one non-redundant peptide set and annotated it once.
# `n_libraries` and `sample_ids` carry the membership through so any tier count
# can be traced back to the libraries that produced it.
#
# INPUTS  (via CLI)
#   --groups        validated group table (group, layout, n_libraries, sample_ids)
#   --results-root  results/{run}/02_annotation
#   --assemblers    comma-separated list actually enabled for this run
# OUTPUT
#   --out-manifest  annotation_manifest.tsv (+ .sha256)
# =============================================================================
# NOTE: no `from __future__ import annotations` — see project gotchas.
# =============================================================================

import argparse
import csv
import hashlib
import sys
from collections import defaultdict
from pathlib import Path


MANIFEST_COLS = [
    "group", "layout", "n_libraries", "sample_ids",
    # R2 ORF prediction (summed over the enabled assemblers)
    "n_assemblers", "n_orfs_total", "n_blast_support_hits",
    # R3 dereplication
    "n_input_orfs", "n_nr_orfs", "n_clusters",
    # R4 ConoSorter
    "n_regex_hits", "n_phmm_hits", "n_regex_and_phmm", "n_candidates",
    # R5 SignalP
    "signalp_route", "n_signalp_positive", "n_after_containment",
    "n_assigned", "n_unassigned",
    # R6 DIAMOND
    "n_blast_hits_uniprot", "n_blast_hits_ncbi", "n_blast_hits_conoserver",
    # R7 labels
    "n_labels_conoserver", "n_labels_ncbi", "n_labels_uniprot",
    # R8 normalization
    "n_normalized_rows", "n_unresolved_labels",
    # R9 tiers
    "n_alta", "n_media", "n_baja", "n_sin_clasificar",
    # provenance
    "classification_xlsx", "notes",
]


def read_kv(p: Path):
    if not p.exists():
        return {}
    out = {}
    with p.open(encoding="utf-8", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for row in rd:
            if len(row) >= 2:
                out[row[0].strip()] = row[1].strip()
    return out


def read_table(p: Path):
    if not p.exists():
        return []
    with p.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def to_int(v, default=0):
    try:
        return int(float(str(v).strip()))
    except (TypeError, ValueError):
        return default


def collect(root: Path, grec: dict, assemblers: list):
    gid = grec["group"]
    base = root / gid
    row = {c: "" for c in MANIFEST_COLS}
    row.update({
        "group": gid,
        "layout": grec.get("layout", ""),
        "n_libraries": grec.get("n_libraries", ""),
        "sample_ids": grec.get("sample_ids", ""),
    })
    notes = []

    # R2 — one predict dir per assembler
    n_asm = n_orfs = n_support = 0
    for asm in assemblers:
        s = read_kv(base / asm / "predict" / "predict.summary.tsv")
        if s:
            n_asm += 1
            n_orfs += to_int(s.get("n_orfs"))
            n_support += to_int(s.get("n_blast_hits"))
    row["n_assemblers"] = n_asm
    row["n_orfs_total"] = n_orfs
    row["n_blast_support_hits"] = n_support
    if n_asm < len(assemblers):
        notes.append(f"predict_missing_for_{len(assemblers) - n_asm}_assembler(s)")

    # R3
    d = read_kv(base / "dedup" / "dedup.summary.tsv")
    row["n_input_orfs"] = to_int(d.get("n_input_orfs"))
    row["n_nr_orfs"] = to_int(d.get("n_nr_orfs"))
    row["n_clusters"] = to_int(d.get("n_clusters"))

    # R4
    m = read_kv(base / "conosorter" / "merge.audit.tsv")
    row["n_regex_hits"] = to_int(m.get("n_regex_hits"))
    row["n_phmm_hits"] = to_int(m.get("n_phmm_hits"))
    row["n_regex_and_phmm"] = to_int(m.get("n_regex_and_phmm"))
    row["n_candidates"] = to_int(m.get("n_unique_seqs"))

    # R5
    s = read_kv(base / "signalp" / "signalp_prepare.audit.tsv")
    row["signalp_route"] = s.get("signalp_route", "")
    row["n_signalp_positive"] = to_int(s.get("n_signalp_positive"))
    row["n_after_containment"] = to_int(s.get("n_after_containment"))
    row["n_assigned"] = to_int(s.get("n_assigned"))
    row["n_unassigned"] = to_int(s.get("n_unassigned"))

    # R6
    per_db = defaultdict(int)
    for r in read_table(base / "annotate" / "blast" / "diamond.summary.tsv"):
        per_db[r.get("db", "")] += to_int(r.get("n_hits"))
    for db in ("uniprot", "ncbi", "conoserver"):
        row[f"n_blast_hits_{db}"] = per_db.get(db, 0)

    # R7
    a = read_kv(base / "annotate" / "annotate_audit.tsv")
    for db in ("conoserver", "ncbi", "uniprot"):
        row[f"n_labels_{db}"] = (to_int(a.get(f"filled_positivo_{db}"))
                                 + to_int(a.get(f"filled_negativo_{db}")))

    # R8
    n = read_kv(base / "classify" / "normalize_audit.tsv")
    row["n_normalized_rows"] = to_int(n.get("rows_total"))
    row["n_unresolved_labels"] = to_int(n.get("unresolved_rows"))

    # R9
    tiers = {"ALTA": 0, "MEDIA": 0, "BAJA": 0, "SIN_CLASIFICAR": 0}
    for r in read_table(base / "classify" / "tier_summary.tsv"):
        if r.get("part") == "all":
            tiers[r.get("tier", "")] = to_int(r.get("n"))
    row["n_alta"] = tiers["ALTA"]
    row["n_media"] = tiers["MEDIA"]
    row["n_baja"] = tiers["BAJA"]
    row["n_sin_clasificar"] = tiers["SIN_CLASIFICAR"]

    xlsx = base / "classify" / "classification.xlsx"
    row["classification_xlsx"] = str(xlsx) if xlsx.exists() else ""
    if not xlsx.exists():
        notes.append("classification_xlsx_missing")

    # A real library can legitimately yield no secreted candidates (wrong
    # tissue, failed prep). Flag it rather than failing the run.
    if row["n_candidates"] and not row["n_signalp_positive"]:
        notes.append("zero_signalp_positive")
    if not row["n_candidates"]:
        notes.append("zero_conosorter_candidates")

    row["notes"] = ";".join(dict.fromkeys(notes))
    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--groups",       required=True, type=Path)
    ap.add_argument("--results-root", required=True, type=Path)
    ap.add_argument("--assemblers",   required=True)
    ap.add_argument("--out-manifest", required=True, type=Path)
    args = ap.parse_args()

    assemblers = [a.strip() for a in args.assemblers.split(",") if a.strip()]
    groups = read_table(args.groups)
    if not groups:
        print(f"[annotation_manifest] ERROR: no groups in {args.groups}",
              file=sys.stderr)
        return 1

    rows = [collect(args.results_root, g, assemblers) for g in groups]

    args.out_manifest.parent.mkdir(parents=True, exist_ok=True)
    with args.out_manifest.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    digest = hashlib.sha256(args.out_manifest.read_bytes()).hexdigest()
    Path(str(args.out_manifest) + ".sha256").write_text(
        f"{digest}  {args.out_manifest.name}\n")

    tot = {k: sum(to_int(r[k]) for r in rows)
           for k in ("n_alta", "n_media", "n_baja", "n_sin_clasificar")}
    flagged = [r["group"] for r in rows if r["notes"]]
    print(f"[annotation_manifest] {len(rows)} group(s); tiers={tot}",
          file=sys.stderr)
    if flagged:
        print(f"[annotation_manifest] flagged: {', '.join(flagged)}",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
