#!/usr/bin/env python3
# =============================================================================
# 07_build_assembly_manifest.py  -  Stage 01 output contract
# =============================================================================
# One row per (group, assembler): what was assembled, how it was filtered, and
# where the file Stage 02 must read lives.
#
# The row key is the GROUP, not the library, because the group is this stage's
# unit of work: Stage 00 concatenated the group's cleaned libraries and Stage
# 01 assembled that pooled read set once. `n_libraries` and `sample_ids`
# carry the membership through, so a reader can trace any assembly back to the
# libraries that produced it without opening a second file.
#
# `filter_method` is the column that matters most downstream. It is
# `transrate` for paired-end libraries and `none` for single-end ones (and for
# any paired sample where TransRate failed and the raw contigs were passed
# through). Contig counts are NOT comparable across those two groups: the
# unfiltered sets retain low-support contigs that TransRate would have
# removed. Any cross-sample analysis has to stratify on it.
#
# INPUTS  (via CLI)
#   --groups        validated group table (group, layout, n_libraries, sample_ids)
#   --assembly-dir  results/{run}/01_assembly
#   --assemblers    comma-separated, e.g. trinity,rnaspades,rnabloom
# OUTPUT
#   --out-manifest  assembly_manifest.tsv (+ .sha256)
# =============================================================================
# NOTE: no `from __future__ import annotations` — see project gotchas.
# =============================================================================

import argparse
import csv
import hashlib
import sys
from pathlib import Path


MANIFEST_COLS = [
    "group", "assembler", "layout", "n_libraries", "sample_ids",
    "assembler_version", "container_image",
    "threads", "max_memory_gb", "assembler_elapsed_s",
    "n_contigs", "n_bases",
    "filter_method", "transrate_score",
    "n_contigs_filtered", "n_bases_filtered", "pct_contigs_retained",
    "transrate_elapsed_s",
    "contigs_fasta", "filtered_contigs_fasta",
    "notes",
]


def read_kv(path: Path):
    if not path.exists():
        return {}
    out = {}
    with path.open(encoding="utf-8", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for row in rd:
            if len(row) >= 2:
                out[row[0].strip()] = row[1].strip()
    return out


def to_int(v, default=0):
    try:
        return int(float(str(v).strip()))
    except (TypeError, ValueError):
        return default


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--groups",       required=True, type=Path)
    ap.add_argument("--assembly-dir", required=True, type=Path)
    ap.add_argument("--assemblers",   required=True)
    ap.add_argument("--min-contigs",  type=int, default=1)
    ap.add_argument("--out-manifest", required=True, type=Path)
    args = ap.parse_args()

    assemblers = [a.strip() for a in args.assemblers.split(",") if a.strip()]

    with args.groups.open(encoding="utf-8", newline="") as fh:
        groups = list(csv.DictReader(fh, delimiter="\t"))
    if not groups:
        print(f"[assembly_manifest] ERROR: no groups in {args.groups}",
              file=sys.stderr)
        return 1

    rows = []
    for grec in groups:
        gid = grec["group"]
        layout = grec.get("layout", "")
        n_libs = grec.get("n_libraries", "")
        sids = grec.get("sample_ids", "")
        for asm in assemblers:
            base = args.assembly_dir / gid / asm
            a = read_kv(base / "run.summary.tsv")
            t = read_kv(base / "transrate.summary.tsv")
            notes = []

            if not a:
                notes.append("assembler_summary_missing")
            if not t:
                notes.append("transrate_summary_missing")

            n_contigs = to_int(a.get("n_contigs"))
            n_filt = to_int(t.get("n_contigs_filtered"))

            if n_contigs < args.min_contigs:
                notes.append(f"below_min_contigs_{args.min_contigs}")
            if n_contigs and n_filt == 0:
                notes.append("all_contigs_filtered_out")
            for extra in (t.get("notes") or "").split(";"):
                if extra:
                    notes.append(extra)

            rows.append({
                "group": gid,
                "assembler": asm,
                "layout": layout,
                "n_libraries": n_libs,
                "sample_ids": sids,
                "assembler_version": a.get("assembler_version", ""),
                "container_image": a.get("container_image", ""),
                "threads": a.get("threads", ""),
                "max_memory_gb": a.get("max_memory_gb", ""),
                "assembler_elapsed_s": a.get("elapsed_seconds", ""),
                "n_contigs": n_contigs,
                "n_bases": to_int(a.get("n_bases")),
                "filter_method": t.get("filter_method", ""),
                "transrate_score": t.get("transrate_score", ""),
                "n_contigs_filtered": n_filt,
                "n_bases_filtered": to_int(t.get("n_bases_filtered")),
                "pct_contigs_retained": t.get("pct_contigs_retained", ""),
                "transrate_elapsed_s": t.get("elapsed_seconds", ""),
                "contigs_fasta": a.get("contigs_fasta", ""),
                "filtered_contigs_fasta": t.get("filtered_fasta", ""),
                "notes": ";".join(dict.fromkeys(notes)),
            })

    args.out_manifest.parent.mkdir(parents=True, exist_ok=True)
    with args.out_manifest.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    digest = hashlib.sha256(args.out_manifest.read_bytes()).hexdigest()
    Path(str(args.out_manifest) + ".sha256").write_text(
        f"{digest}  {args.out_manifest.name}\n")

    by_method = {}
    for r in rows:
        by_method[r["filter_method"]] = by_method.get(r["filter_method"], 0) + 1
    flagged = [f"{r['group']}/{r['assembler']}" for r in rows if r["notes"]]

    print(f"[assembly_manifest] {len(rows)} rows "
          f"({len(groups)} groups x {len(assemblers)} assemblers); "
          f"filter_method={by_method}", file=sys.stderr)
    if flagged:
        print(f"[assembly_manifest] flagged: {', '.join(flagged[:15])}"
              + (" ..." if len(flagged) > 15 else ""), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
