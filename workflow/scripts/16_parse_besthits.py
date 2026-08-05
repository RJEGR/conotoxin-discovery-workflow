#!/usr/bin/env python3
# =============================================================================
# 29_stage7_parse_diamond_besthits.py  -  Stage 7 R6b: outfmt6 -> TSV +
#                                          best-hit selection + coverage audit.
# =============================================================================
# Port of legacy `7_diamond_outfmt6_to_tsv_pipeline.py`, adapted to the new
# per-fold directory layout and to emit a single consolidated best-hits table
# that R7 consumes (the legacy version left best-hit selection to scripts
# 8/9/10, duplicating the logic three times).
#
# Best-hit rule (verbatim from <context> annotation step 3): "the best-supported
# blastp hits were selected on the basis of the highest bit score and lowest
# e-value".
#
# Per-DB hit retention differs, mirroring the legacy behaviour:
#   * conoserver -> ALL hits retained (--max-target-seqs 3); labels are later
#     collapsed comma-joined by R7. Rationale: ConoServer entries are
#     fine-grained conopeptide records, and several may map to the same
#     superfamily, so keeping them all lets R7 report agreement.
#   * uniprot / ncbi -> single best hit per query (--max-target-seqs 2 gives
#     the parser a runner-up to audit against, but only the best is used).
#
# Coverage audit (legacy "flagged_alt_more_informative"): flags queries where
# the top-bitscore hit has poor query coverage but an alternative hit with a
# comparable bitscore covers substantially more of the query. These are not
# corrected automatically -- they are reported so a curator can inspect them.
#
# INPUTS (via CLI):
#   --blast-dir    dir holding {part}__{db}.outfmt6  (6 files)
#   --min-qcov     float, default 30.0
#   --cov-delta    float, default 20.0
#   --bitscore-ratio float, default 0.80
#
# OUTPUTS (via CLI):
#   --out-besthits  besthits.tsv    (qseqid, part, db, sseqid, accession,
#                                    bitscore, evalue, qcovhsp, n_hits, all_sseqids)
#   --out-tsv-dir   dir for the 6 header-ed .tsv conversions
#   --out-audit     blast_audit.tsv (per-file counts + flagged rows summary)
#   --out-flagged   blast_flagged_alt_more_informative.tsv
# =============================================================================
# NOTE: no `from __future__ import annotations` — see gotcha #1.
# =============================================================================

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


COLS = [
    "qseqid", "qlen", "sseqid", "slen", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "scovhsp",
]
DB_TAGS = ("uniprot", "ncbi", "conoserver")
PARTS = ("positivo", "negativo")

# DBs for which every hit is retained (labels collapsed downstream by R7).
MULTI_HIT_DBS = {"conoserver"}


def safe_float(x, default=float("nan")):
    try:
        return float(str(x).strip())
    except (TypeError, ValueError):
        return default


def normalize_accession(sseqid: str, db: str) -> str:
    """Extract the DB-native accession from a DIAMOND sseqid.

    ConoServer : `P00001|SI|Conus striatus|...`      -> `P00001`  (field 0)
    UniProt    : `sp|A0A068B6Q6|CA18_CONBE`          -> `A0A068B6Q6` (field 1)
    NCBI       : `sp|P0DRL4.1|CAIA_CONNO` or bare    -> field 1 if piped, else raw
    """
    s = str(sseqid or "").strip()
    if not s:
        return ""
    if db == "conoserver":
        return s.split("|", 1)[0].strip()
    if "|" in s:
        parts = s.split("|")
        if len(parts) >= 2 and parts[1].strip():
            return parts[1].strip()
    return s.split()[0]


def read_outfmt6(path: Path):
    """Return list of dict rows. Tolerates empty/missing files."""
    rows = []
    if not path.exists() or path.stat().st_size == 0:
        return rows
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) != len(COLS):
                parts = line.split()
            if len(parts) != len(COLS):
                # Malformed row: skip but keep going (legacy passed it through;
                # here we drop it because downstream expects a strict schema).
                continue
            rows.append(dict(zip(COLS, parts)))
    return rows


def write_tsv(rows: list, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in COLS) + "\n")


def audit_coverage(rows: list, min_qcov: float, cov_delta: float,
                   bitscore_ratio: float, part: str, db: str):
    """Legacy `flagged_alt_more_informative` logic. Returns (flagged, stats)."""
    hits_by_q = defaultdict(list)
    for r in rows:
        hits_by_q[r["qseqid"]].append(r)

    flagged = []
    n_qcov_ge_min = 0

    for q, hits in hits_by_q.items():
        enriched = [(safe_float(h["bitscore"], float("-inf")),
                     safe_float(h["qcovhsp"]), h) for h in hits]
        # top1 = highest bitscore, tie-broken by qcov
        top1_bits, top1_qcov, top1 = sorted(
            enriched, key=lambda x: (x[0], x[1] if x[1] == x[1] else -1), reverse=True)[0]
        # bestcov = highest qcov, tie-broken by bitscore
        bc_bits, bc_qcov, bc = sorted(
            enriched, key=lambda x: (x[1] if x[1] == x[1] else -1, x[0]), reverse=True)[0]

        if any(e[1] == e[1] and e[1] >= min_qcov for e in enriched):
            n_qcov_ge_min += 1

        reason = None
        if (top1_qcov != top1_qcov) or (top1_qcov < min_qcov):
            if (bc_qcov == bc_qcov) and (bc_qcov >= min_qcov):
                if (top1_bits == top1_bits) and (bc_bits >= bitscore_ratio * top1_bits):
                    reason = (f"top1_qcov<{min_qcov} but alt_qcov>={min_qcov} "
                              f"and bitscore_ratio>={bitscore_ratio}")
                else:
                    reason = f"top1_qcov<{min_qcov} but alt_qcov>={min_qcov}"
        if reason is None and (top1_qcov == top1_qcov) and (bc_qcov == bc_qcov):
            if (bc_qcov >= top1_qcov + cov_delta) and (bc_bits >= bitscore_ratio * top1_bits):
                reason = (f"alt_qcov>=top1_qcov+{cov_delta} "
                          f"and bitscore_ratio>={bitscore_ratio}")

        if reason and (bc["sseqid"] != top1["sseqid"] or bc_qcov != top1_qcov):
            flagged.append({
                "part": part, "db": db, "qseqid": q,
                "top1_sseqid": top1["sseqid"], "top1_bitscore": top1_bits,
                "top1_evalue": top1["evalue"], "top1_qcovhsp": top1["qcovhsp"],
                "bestcov_sseqid": bc["sseqid"], "bestcov_bitscore": bc_bits,
                "bestcov_evalue": bc["evalue"], "bestcov_qcovhsp": bc["qcovhsp"],
                "reason": reason,
            })

    stats = {
        "part": part, "db": db,
        "n_rows": len(rows),
        "n_queries_with_hits": len(hits_by_q),
        "n_queries_qcov_ge_min": n_qcov_ge_min,
        "n_flagged": len(flagged),
        "pct_flagged": round(len(flagged) / len(hits_by_q) * 100, 2) if hits_by_q else 0.0,
    }
    return flagged, stats


def select_besthits(rows: list, part: str, db: str):
    """Best hit per query: highest bitscore, then lowest e-value, then input order.

    For MULTI_HIT_DBS we additionally record every sseqid (semicolon-joined)
    so R7 can collapse their labels.
    """
    hits_by_q = defaultdict(list)
    for i, r in enumerate(rows):
        hits_by_q[r["qseqid"]].append((i, r))

    out = []
    for q, indexed in hits_by_q.items():
        ordered = sorted(
            indexed,
            key=lambda ir: (-safe_float(ir[1]["bitscore"], float("-inf")),
                            safe_float(ir[1]["evalue"], float("inf")),
                            ir[0]),
        )
        best = ordered[0][1]
        all_sseqids = ";".join(r["sseqid"] for _i, r in ordered)
        out.append({
            "qseqid": q,
            "part": part,
            "db": db,
            "sseqid": best["sseqid"],
            "accession": normalize_accession(best["sseqid"], db),
            "bitscore": best["bitscore"],
            "evalue": best["evalue"],
            "qcovhsp": best["qcovhsp"],
            "pident": best["pident"],
            "n_hits": len(ordered),
            "all_sseqids": all_sseqids if db in MULTI_HIT_DBS else best["sseqid"],
        })
    return out


BESTHIT_COLS = ["qseqid", "part", "db", "sseqid", "accession", "bitscore",
                "evalue", "qcovhsp", "pident", "n_hits", "all_sseqids"]

AUDIT_COLS = ["part", "db", "n_rows", "n_queries_with_hits",
              "n_queries_qcov_ge_min", "n_flagged", "pct_flagged"]

FLAGGED_COLS = ["part", "db", "qseqid",
                "top1_sseqid", "top1_bitscore", "top1_evalue", "top1_qcovhsp",
                "bestcov_sseqid", "bestcov_bitscore", "bestcov_evalue",
                "bestcov_qcovhsp", "reason"]


def write_table(rows: list, cols: list, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blast-dir",       required=True, type=Path)
    ap.add_argument("--out-besthits",    required=True, type=Path)
    ap.add_argument("--out-tsv-dir",     required=True, type=Path)
    ap.add_argument("--out-audit",       required=True, type=Path)
    ap.add_argument("--out-flagged",     required=True, type=Path)
    ap.add_argument("--min-qcov",        type=float, default=30.0)
    ap.add_argument("--cov-delta",       type=float, default=20.0)
    ap.add_argument("--bitscore-ratio",  type=float, default=0.80)
    args = ap.parse_args()

    all_besthits = []
    all_audit = []
    all_flagged = []

    for part in PARTS:
        for db in DB_TAGS:
            src = args.blast_dir / f"{part}__{db}.outfmt6"
            rows = read_outfmt6(src)
            write_tsv(rows, args.out_tsv_dir / f"{part}__{db}.tsv")
            flagged, stats = audit_coverage(rows, args.min_qcov, args.cov_delta,
                                            args.bitscore_ratio, part, db)
            all_flagged.extend(flagged)
            all_audit.append(stats)
            all_besthits.extend(select_besthits(rows, part, db))

    write_table(all_besthits, BESTHIT_COLS, args.out_besthits)
    write_table(all_audit,    AUDIT_COLS,   args.out_audit)
    write_table(all_flagged,  FLAGGED_COLS, args.out_flagged)

    n_q = len({(r["qseqid"], r["part"]) for r in all_besthits})
    print(f"[parse_besthits] besthits={len(all_besthits)} rows over "
          f"{n_q} (query,part) pairs; flagged={len(all_flagged)}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
