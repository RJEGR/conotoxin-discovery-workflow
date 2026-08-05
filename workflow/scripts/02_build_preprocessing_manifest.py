#!/usr/bin/env python3
# =============================================================================
# 02_build_preprocessing_manifest.py  -  Stage 00 output contract
# =============================================================================
# Collapses the per-sample fastp summaries into one flat TSV: the contract
# Stage 01 consumes and the record of what preprocessing actually did.
#
# The manifest is emitted in BOTH modes:
#   * fastp enabled  -> one row per sample with the full before/after metrics
#     and `reads_source = fastp`.
#   * fastp disabled -> one row per sample pointing at the raw FASTQs, with
#     `reads_source = raw` and the metric columns left blank.
#
# Emitting it either way is deliberate: downstream stages and any archived run
# then have a uniform provenance record, and flipping the toggle changes a
# value in the table rather than the shape of the results tree.
#
# INPUTS  (via CLI)
#   --samples        validated sample sheet (sample_id, fastq_1, fastq_2, layout)
#   --summary-dir    directory holding {sample}/{sample}_fastp.summary.tsv
#   --fastp-enabled  true | false
#
# OUTPUT
#   --out-manifest   preprocessing_manifest.tsv  (+ .sha256 sidecar)
# =============================================================================
# NOTE: no `from __future__ import annotations` — see project gotchas.
# =============================================================================

import argparse
import csv
import hashlib
import sys
from pathlib import Path


MANIFEST_COLS = [
    "sample_id", "group", "layout", "reads_source",
    "group_n_libraries", "group_reads_r1", "group_reads_r2",
    "reads_r1", "reads_r2",
    "raw_fastq_1", "raw_fastq_2",
    "reads_before", "reads_after", "pct_reads_passed",
    "bases_before", "bases_after",
    "q20_rate_after", "q30_rate_after", "gc_content_after",
    "read1_mean_length_after", "read2_mean_length_after",
    "duplication_rate", "reads_low_quality", "reads_too_many_N",
    "reads_too_short",
    "fastp_version", "elapsed_seconds", "notes",
]

# Columns copied verbatim from the fastp summary when it is present.
_FROM_SUMMARY = [
    "reads_before", "reads_after", "pct_reads_passed",
    "bases_before", "bases_after",
    "q20_rate_after", "q30_rate_after", "gc_content_after",
    "read1_mean_length_after", "read2_mean_length_after",
    "duplication_rate", "reads_low_quality", "reads_too_many_N",
    "reads_too_short", "fastp_version", "elapsed_seconds",
]


def read_kv(path: Path):
    """Two-column metric/value TSV -> dict. Missing file -> {}."""
    if not path.exists():
        return {}
    out = {}
    with path.open(encoding="utf-8", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)                       # header
        for row in rd:
            if len(row) >= 2:
                out[row[0].strip()] = row[1].strip()
    return out


def to_float(v, default=None):
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return default


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples",       required=True, type=Path)
    ap.add_argument("--summary-dir",   required=True, type=Path)
    ap.add_argument("--group-dir",     required=False, type=Path, default=None)
    ap.add_argument("--fastp-enabled", default="true")
    ap.add_argument("--out-manifest",  required=True, type=Path)
    args = ap.parse_args()

    fastp_on = str(args.fastp_enabled).strip().lower() == "true"

    with args.samples.open(encoding="utf-8", newline="") as fh:
        samples = list(csv.DictReader(fh, delimiter="\t"))
    if not samples:
        print(f"[preprocessing_manifest] ERROR: no samples in {args.samples}",
              file=sys.stderr)
        return 1

    # Pooling summaries, read once and shared by every member of a group.
    group_info = {}
    if args.group_dir:
        for s in samples:
            g = (s.get("group") or s["sample_id"]).strip()
            if g not in group_info:
                group_info[g] = read_kv(
                    args.group_dir / g / f"{g}_concat.summary.tsv")

    rows = []
    for s in samples:
        sid = s["sample_id"]
        grp = (s.get("group") or sid).strip()
        layout = s.get("layout") or ("PE" if s.get("fastq_2") else "SE")
        notes = []

        row = {c: "" for c in MANIFEST_COLS}
        row.update({
            "sample_id": sid,
            "group": grp,
            "layout": layout,
            "raw_fastq_1": s.get("fastq_1", ""),
            "raw_fastq_2": s.get("fastq_2", ""),
        })

        # Group-level columns are repeated on every member row so the manifest
        # stays a flat table: one row per library, joinable on `group` without
        # a second file.
        gi = group_info.get(grp, {})
        row["group_n_libraries"] = gi.get("n_libraries", "")
        row["group_reads_r1"] = gi.get("concat_r1", "")
        row["group_reads_r2"] = gi.get("concat_r2", "")
        if args.group_dir and not gi:
            notes.append("group_concat_summary_missing")

        if fastp_on:
            summ = read_kv(args.summary_dir / sid / f"{sid}_fastp.summary.tsv")
            if not summ:
                notes.append("fastp_summary_missing")
            row["reads_source"] = "fastp"
            row["reads_r1"] = summ.get("clean_r1", "")
            row["reads_r2"] = summ.get("clean_r2", "")
            for c in _FROM_SUMMARY:
                row[c] = summ.get(c, "")

            # QC flags. These are advisory: a low retention rate can be
            # legitimate for a degraded library, but it should be visible in
            # the contract rather than buried in an HTML report.
            pct = to_float(summ.get("pct_reads_passed"))
            if pct is not None and pct < 50.0:
                notes.append(f"low_retention_{pct:.1f}pct")
            after = to_float(summ.get("reads_after"))
            if after is not None and after == 0:
                notes.append("zero_reads_after_filtering")

            # Cross-check fastp's own count against the file actually written.
            # For PE, fastp reports total reads across BOTH mates while seqkit
            # counts records in R1 only, so the expected ratio is 2:1.
            sk = to_float(summ.get("seqkit_reads_out_r1"))
            if after is not None and sk is not None and sk > 0:
                expected = after / 2 if layout == "PE" else after
                if abs(expected - sk) > max(1.0, 0.001 * expected):
                    notes.append(
                        f"count_mismatch_json{int(expected)}_vs_file{int(sk)}")
        else:
            row["reads_source"] = "raw"
            row["reads_r1"] = s.get("fastq_1", "")
            row["reads_r2"] = s.get("fastq_2", "")
            notes.append("fastp_disabled")

        row["notes"] = ";".join(notes)
        rows.append(row)

    args.out_manifest.parent.mkdir(parents=True, exist_ok=True)
    with args.out_manifest.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    digest = hashlib.sha256(args.out_manifest.read_bytes()).hexdigest()
    Path(str(args.out_manifest) + ".sha256").write_text(
        f"{digest}  {args.out_manifest.name}\n")

    n_pe = sum(1 for r in rows if r["layout"] == "PE")
    n_se = len(rows) - n_pe
    n_groups = len({r["group"] for r in rows})
    n_pooled = len({r["group"] for r in rows
                    if to_float(r["group_n_libraries"], 1) > 1})
    flagged = [r["sample_id"] for r in rows
               if r["notes"] and r["notes"] != "fastp_disabled"]
    print(f"[preprocessing_manifest] {len(rows)} libraries in {n_groups} "
          f"group(s) ({n_pooled} pooling >1 library); {n_pe} PE, {n_se} SE; "
          f"source={'fastp' if fastp_on else 'raw'}", file=sys.stderr)
    if flagged:
        print(f"[preprocessing_manifest] flagged: {', '.join(flagged)}",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
