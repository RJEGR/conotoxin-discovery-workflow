#!/usr/bin/env python3
# =============================================================================
# 00_validate_samples.py  -  validate the input sample sheet
# =============================================================================
# Runs before anything else. A malformed sample sheet is the single most common
# way to waste hours of cluster time, so every check that can be made cheaply
# up front is made here, and the failure message names the offending row.
#
# CHECKS
#   1. File exists, is a TSV, and carries the required header
#      (`sample_id`, `fastq_1`; `group` and `fastq_2` optional).
#   2. `sample_id` is present, unique, and safe as a directory name and
#      Snakemake wildcard: [A-Za-z0-9._-] only.
#   3. `fastq_1` is a non-empty ABSOLUTE path that exists and is readable.
#   4. `fastq_2`, if present and non-empty, gets the same treatment, and must
#      differ from `fastq_1`.
#   5. Layout is derived: PE when fastq_2 is populated, SE otherwise.
#   6. FASTQ extension sanity -- a WARNING, not a failure.
#   7. Files are non-empty (size > 0).
#   8. `group` obeys the same character rules as `sample_id` and defaults to
#      the sample_id when blank or absent.
#   9. Every library in a group shares one layout. Mixing SE and PE inside a
#      group is a hard failure: concatenating unpaired with paired reads makes
#      an incoherent read set that assemblers either reject or silently
#      mis-assemble.
#  10. A group name may not collide with the sample_id of a library that
#      belongs to a DIFFERENT group -- two assembly units would otherwise
#      resolve to the same output directory.
#
# GROUPS: THE ASSEMBLY UNIT
#   Stage 00 (fastp) runs per LIBRARY. Adapter content and quality differ
#   between lanes and sequencing runs, so each library is cleaned with its own
#   detection and keeps its own QC report — that is also how a bad replicate
#   gets spotted.
#   Stages 01 and 02 run per GROUP. The group's cleaned reads are concatenated
#   and assembled once. Assembling replicates separately fragments transcripts
#   that the pooled depth recovers intact, which matters especially for
#   conotoxins: precursors are short and expression is highly uneven, so a
#   replicate on its own often yields partial ORFs that fail the signal-peptide
#   test downstream.
#   With no `group` column (or an empty value) each library is its own group,
#   reproducing one-assembly-per-library behaviour.
#
# OUTPUT
#   samples_validated.tsv : sample_id, group, fastq_1, fastq_2, layout, n_files
#   groups_validated.tsv  : group, layout, n_libraries, sample_ids
#   validation_report.txt : human-readable per-check log
#
# Exits non-zero on any hard failure. This script is intentionally the ONE
# place allowed to hard-fail on input problems, so downstream rules can assume
# a clean sheet.
# =============================================================================
# NOTE: no `from __future__ import annotations` — see project gotchas.
# =============================================================================

import argparse
import csv
import os
import re
import sys
from pathlib import Path


SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9._-]+$")
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz",
                  ".fastq.bz2", ".fq.bz2")

REQUIRED_COLS = ("sample_id", "fastq_1")
OUT_COLS = ["sample_id", "group", "fastq_1", "fastq_2", "layout", "n_files"]
GROUP_COLS = ["group", "layout", "n_libraries", "sample_ids"]


def looks_like_fastq(path: str) -> bool:
    low = path.lower()
    return any(low.endswith(s) for s in FASTQ_SUFFIXES)


def check_read_file(path: str, label: str, sample: str, findings: list) -> bool:
    """Return True if the file passes every hard check."""
    ok = True
    p = Path(path)

    if not p.is_absolute():
        findings.append((sample, label, "FAIL", path,
                         "path is not absolute; use a full path so the sheet "
                         "resolves identically on submit and compute nodes"))
        return False

    if not p.exists():
        findings.append((sample, label, "FAIL", path, "file does not exist"))
        return False

    if not p.is_file():
        findings.append((sample, label, "FAIL", path, "not a regular file"))
        ok = False
    elif p.stat().st_size == 0:
        findings.append((sample, label, "FAIL", path, "file is empty"))
        ok = False
    elif not os.access(p, os.R_OK):
        findings.append((sample, label, "FAIL", path, "file is not readable"))
        ok = False
    else:
        findings.append((sample, label, "PASS", path,
                         f"{p.stat().st_size:,} bytes"))

    if not looks_like_fastq(path):
        findings.append((sample, label, "WARN", path,
                         "unusual extension for a FASTQ; proceeding anyway"))
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples",    required=True, type=Path)
    ap.add_argument("--out-tsv",    required=True, type=Path)
    ap.add_argument("--out-groups", required=False, type=Path, default=None)
    ap.add_argument("--out-report", required=True, type=Path)
    args = ap.parse_args()

    findings = []
    fails = 0
    rows_out = []

    def add(sample, check, status, value, detail):
        nonlocal fails
        findings.append((sample, check, status, value, detail))
        if status == "FAIL":
            fails += 1

    # ---- 1. sheet itself --------------------------------------------------
    if not args.samples.exists():
        add("-", "samples_file", "FAIL", str(args.samples), "does not exist")
        _emit(args, findings, fails, rows_out, [])
        return 1
    if args.samples.stat().st_size == 0:
        add("-", "samples_file", "FAIL", str(args.samples), "is empty")
        _emit(args, findings, fails, rows_out, [])
        return 1
    add("-", "samples_file", "PASS", str(args.samples),
        f"{args.samples.stat().st_size} bytes")

    with args.samples.open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = [h.strip() for h in (reader.fieldnames or [])]
        missing = [c for c in REQUIRED_COLS if c not in header]
        if missing:
            add("-", "header", "FAIL", ",".join(header),
                f"missing required column(s): {missing}. "
                f"Expected at least {list(REQUIRED_COLS)}. "
                f"Columns must be separated by TAB characters, not spaces — "
                f"check with `cat -A config/samples.tsv` (tabs show as ^I).")
            _emit(args, findings, fails, rows_out, [])
            return 1
        add("-", "header", "PASS", ",".join(header),
            "fastq_2 " + ("present" if "fastq_2" in header else
                          "absent -> all libraries treated as single-end"))
        raw_rows = list(reader)

    if not raw_rows:
        add("-", "rows", "FAIL", "0", "sample sheet has a header but no rows")
        _emit(args, findings, fails, rows_out, [])
        return 1
    add("-", "rows", "PASS", str(len(raw_rows)), "data rows found")

    # ---- 2..7 per-row checks ---------------------------------------------
    seen_ids = {}
    for i, row in enumerate(raw_rows, start=2):     # line 1 is the header
        sid = (row.get("sample_id") or "").strip()
        grp = (row.get("group") or "").strip()
        f1 = (row.get("fastq_1") or "").strip()
        f2 = (row.get("fastq_2") or "").strip()

        if not sid:
            add(f"<line {i}>", "sample_id", "FAIL", "", "sample_id is empty")
            continue
        if not SAMPLE_ID_RE.match(sid):
            add(sid, "sample_id", "FAIL", sid,
                "contains characters outside [A-Za-z0-9._-]; it is used as a "
                "directory name and a Snakemake wildcard")
            continue
        if sid in seen_ids:
            add(sid, "sample_id", "FAIL", sid,
                f"duplicated (first seen on line {seen_ids[sid]})")
            continue
        seen_ids[sid] = i
        add(sid, "sample_id", "PASS", sid, f"line {i}")

        # ---- 8. group (assembly unit) ------------------------------------
        if not grp:
            grp = sid
            add(sid, "group", "PASS", grp,
                "no group given -> library assembles on its own")
        elif not SAMPLE_ID_RE.match(grp):
            add(sid, "group", "FAIL", grp,
                "contains characters outside [A-Za-z0-9._-]; the group is used "
                "as a directory name and a Snakemake wildcard")
            continue
        else:
            add(sid, "group", "PASS", grp, "assembly unit")

        if not f1:
            add(sid, "fastq_1", "FAIL", "", "fastq_1 is empty")
            continue
        ok1 = check_read_file(f1, "fastq_1", sid, findings)
        if not ok1:
            fails += 1

        layout = "SE"
        ok2 = True
        if f2:
            if f2 == f1:
                add(sid, "fastq_2", "FAIL", f2,
                    "fastq_2 is identical to fastq_1")
                ok2 = False
            else:
                ok2 = check_read_file(f2, "fastq_2", sid, findings)
                if not ok2:
                    fails += 1
                layout = "PE"
        else:
            add(sid, "layout", "PASS", "SE",
                "fastq_2 empty -> single-end. NOTE: TransRate filtering is "
                "skipped for SE libraries (see config assembly.transrate)")

        if ok1 and ok2:
            rows_out.append({
                "sample_id": sid,
                "group": grp,
                "fastq_1": f1,
                "fastq_2": f2,
                "layout": layout,
                "n_files": 2 if layout == "PE" else 1,
            })

    # ---- 9/10. group-level consistency ------------------------------------
    groups = {}
    for r in rows_out:
        groups.setdefault(r["group"], []).append(r)

    group_rows = []
    for gname, members in sorted(groups.items()):
        layouts = sorted({m["layout"] for m in members})
        sids = [m["sample_id"] for m in members]

        # 9. one layout per group
        if len(layouts) > 1:
            detail = "; ".join(f"{m['sample_id']}={m['layout']}"
                               for m in members)
            add(gname, "group_layout", "FAIL", "/".join(layouts),
                "group mixes single- and paired-end libraries, which cannot be "
                f"concatenated into one coherent read set ({detail}). Split the "
                "group, or drop the odd libraries.")
            continue

        # 10. a group name may not shadow a foreign sample_id
        if gname in seen_ids and gname not in sids:
            add(gname, "group_name", "FAIL", gname,
                f"group name collides with the sample_id on line "
                f"{seen_ids[gname]}, which belongs to a different group; both "
                "would resolve to the same output directory")
            continue

        add(gname, "group", "PASS", layouts[0],
            f"{len(members)} librar{'y' if len(members) == 1 else 'ies'}: "
            + ", ".join(sids))
        group_rows.append({
            "group": gname,
            "layout": layouts[0],
            "n_libraries": len(members),
            "sample_ids": ",".join(sids),
        })

    # ---- summary ----------------------------------------------------------
    n_pe = sum(1 for r in rows_out if r["layout"] == "PE")
    n_se = sum(1 for r in rows_out if r["layout"] == "SE")
    n_pooled = sum(1 for g in group_rows if g["n_libraries"] > 1)
    add("-", "summary", "PASS",
        f"{len(rows_out)} libraries / {len(group_rows)} groups",
        f"{n_pe} PE, {n_se} SE; {n_pooled} group(s) pool >1 library")

    _emit(args, findings, fails, rows_out, group_rows)

    if fails:
        print(f"[validate_samples] FAILED with {fails} error(s); see "
              f"{args.out_report}", file=sys.stderr)
        return 1
    print(f"[validate_samples] OK: {len(rows_out)} libraries in "
          f"{len(group_rows)} assembly group(s) ({n_pe} PE, {n_se} SE)",
          file=sys.stderr)
    return 0


def _emit(args, findings, fails, rows_out, group_rows):
    args.out_report.parent.mkdir(parents=True, exist_ok=True)
    with args.out_report.open("w", encoding="utf-8") as fh:
        fh.write("# Sample sheet validation\n")
        fh.write(f"# sheet = {args.samples}\n")
        fh.write("# Stage 00 runs per LIBRARY; Stages 01-02 run per GROUP.\n\n")
        w = max((len(f[0]) for f in findings), default=10) + 2
        for sample, check, status, value, detail in findings:
            fh.write(f"[{status:6}] {sample:<{w}} {check:<13} :: {detail}\n")
        fh.write(f"\nTotal FAIL: {fails}\n")

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open("w", encoding="utf-8", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t")
        wr.writeheader()
        wr.writerows(rows_out)

    if args.out_groups:
        args.out_groups.parent.mkdir(parents=True, exist_ok=True)
        with args.out_groups.open("w", encoding="utf-8", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=GROUP_COLS, delimiter="\t")
            wr.writeheader()
            wr.writerows(group_rows)


if __name__ == "__main__":
    sys.exit(main())
