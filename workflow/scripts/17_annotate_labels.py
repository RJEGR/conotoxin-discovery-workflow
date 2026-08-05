#!/usr/bin/env python3
# =============================================================================
# 30_stage7_annotate_all_labels.py  -  Stage 7 R7: resolve annotation labels
#                                       from ConoServer, NCBI and UniProt.
# =============================================================================
# Consolidated replacement for legacy scripts 8, 9 and 10, which each
# implemented the same skeleton (read best-hits -> resolve accession against a
# label source -> write a column into the xlsx) against a different label
# source. Fusing them removes ~370 lines of duplicated best-hit / xlsx /
# audit plumbing while preserving each source's resolution rules verbatim.
#
# LABEL SOURCES (one per DB; each keeps its legacy semantics)
#
#   conoserver  (legacy 8)
#     Source: the ConoServer protein FASTA. Headers are pipe-delimited:
#       >P00001|SI|Conus striatus|Wild type|conotoxin|A superfamily|...
#     Field 6 (1-based, configurable) carries the superfamily. Because
#     ConoServer runs with --max-target-seqs 3, ALL hits are resolved and
#     their distinct labels comma-joined (legacy `build_qseqid_to_collapsed_values`).
#
#   ncbi        (legacy 9)
#     Source: the NCBI feature table (`base_datos_anotacion_ncbi.txt`), an
#     ASN-tabular file of `>Feature gb|ACC|` blocks. The superfamily is mined
#     from `product` / `note` / `region` lines using a PRIORITY-ORDERED regex
#     table, and the captured token must be in a curated WHITELIST of known
#     superfamily symbols. Only the best hit's accession is resolved.
#
#   uniprot     (legacy 10)
#     Source: `base_datos_anotacion_uniprot.xlsx` (cols: Entry, Protein names,
#     Organism). Accession is taken from `sp|ACC|NAME` -> ACC. Resolution is
#     exact-first, then a canonical fallback for isoform IDs (`C0KYC5-2` ->
#     `C0KYC5`) so a missing isoform record still yields the canonical label.
#     Only the best hit's accession is resolved.
#
# INPUTS (via CLI):
#   --besthits            besthits.tsv                    (from R6b)
#   --pos-xlsx/--neg-xlsx signalp/{positivo,negativo}_superfamilia.xlsx
#   --conoserver-fasta    ConoServer protein FASTA
#   --conoserver-field    int, 1-based header field (default 6)
#   --ncbi-feature-table  base_datos_anotacion_ncbi.txt
#   --uniprot-xlsx        base_datos_anotacion_uniprot.xlsx
#   --uniprot-sheet / --uniprot-entry-col / --uniprot-name-col
#   --key-col             xlsx column holding the query id (default seq_id)
#
# OUTPUTS (via CLI):
#   --out-pos-xlsx/--out-neg-xlsx  annotated workbooks (+3 columns)
#   --out-labels                   labels.tsv  (tidy: seq_id, part, conoserver,
#                                               ncbi, uniprot)
#   --out-audit                    annotate_audit.tsv
# =============================================================================
# NOTE: no `from __future__ import annotations` — see gotcha #1.
# =============================================================================

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import openpyxl


DB_TAGS = ("conoserver", "ncbi", "uniprot")

# ---------------------------------------------------------------------------
# NCBI superfamily mining (verbatim from legacy script 9)
# ---------------------------------------------------------------------------
NCBI_WHITELIST = {
    "A", "B1", "B2", "B3", "C", "D", "E", "F", "G", "H", "I1", "I2", "I3", "J",
    "K", "L", "M", "N", "O1", "O2", "O3", "P", "Q", "R", "S", "T", "U", "V",
    "W", "Y", "Z", "conkunitzin",
}

# (source_field, priority, regex). Lower priority number wins.
NCBI_PATTERNS = [
    ("product", 1, re.compile(r"\bproduct\b.*\bconotoxin\s+superfamily\s+([^\s;,\)\]\}]+)", re.IGNORECASE)),
    ("product", 2, re.compile(r"\bproduct\b.*\bconotoxin\s+precursor\s+superfamily\s+([^\s;,\)\]\}]+)", re.IGNORECASE)),
    ("product", 3, re.compile(r"\bproduct\b.*\bconotoxin\s+precursor\s+([^\s;,\)\]\}]+)", re.IGNORECASE)),
    ("note",    4, re.compile(r"\bnote\b.*\bconotoxin\s+([A-Za-z0-9]+)\s*-\s*superfamily\b", re.IGNORECASE)),
    ("note",    5, re.compile(r"\bnote\b.*\b([A-Za-z0-9]+)\s*-\s*superfamily\b", re.IGNORECASE)),
    ("note",    6, re.compile(r"\bnote\b.*\b([A-Za-z0-9]+)\s+superfamily\b", re.IGNORECASE)),
    ("region",  7, re.compile(r"\bregion\b.*\bconotoxin_([A-Za-z0-9]+)\b", re.IGNORECASE)),
]

ISOFORM_SUFFIX_RE = re.compile(r"^(?P<canonical>.+)-(?P<isoform>[1-9][0-9]*)$")


def _normalize_sf_token(token: str):
    t = str(token or "").strip().strip("\"' \t")
    if not t:
        return "", False
    if t.lower() == "conkunitzin":
        return "conkunitzin", True
    t_up = t.upper()
    if t_up in NCBI_WHITELIST:
        return t_up, True
    return t, False


def _extract_best_sf_from_block(block_lines):
    best = None
    for line in block_lines:
        line_str = line.rstrip("\n")
        for source_field, prio, rx in NCBI_PATTERNS:
            m = rx.search(line_str)
            if not m:
                continue
            sf_norm, in_wl = _normalize_sf_token(m.group(1))
            # Legacy behaviour: a token must be in the whitelist to be accepted.
            if not in_wl:
                continue
            cand = (prio, sf_norm, source_field)
            if best is None or cand[0] < best[0]:
                best = cand
    if best is None:
        return ""
    return best[1]


# ---------------------------------------------------------------------------
# Label-source loaders
# ---------------------------------------------------------------------------
def load_conoserver_labels(fasta_path: Path, field_1based: int, sep: str = "|"):
    """accession -> superfamily string, from pipe-delimited FASTA headers."""
    idx = {}
    dup = set()
    with fasta_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = line[1:].rstrip("\n").split(sep)
            accession = parts[0].strip()
            value = parts[field_1based - 1].strip() if len(parts) >= field_1based else ""
            if accession in idx and idx[accession] != value:
                dup.add(accession)
            idx[accession] = value
    return idx, dup


def load_ncbi_labels(feature_table: Path, wanted_accessions: set):
    """accession -> superfamily, streaming only the blocks we care about."""
    ann = {}
    current_acc = None
    want = False
    block = []

    def finalize():
        nonlocal current_acc, want, block
        if current_acc and want:
            ann[current_acc] = _extract_best_sf_from_block(block)
        current_acc, want, block = None, False, []

    with feature_table.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">Feature "):
                finalize()
                token = line.strip().split(" ", 1)[1].strip()
                # `gb|YFQ72895.1|` -> `YFQ72895.1`
                acc = token.split("|")[1].strip() if "|" in token else token
                current_acc = acc
                want = acc in wanted_accessions
                block = [line] if want else []
            elif want:
                block.append(line)
        finalize()
    return ann


def load_uniprot_labels(xlsx_path: Path, sheet: str, entry_col: str, name_col: str):
    """Entry -> Protein names, read with openpyxl (read-only, fast)."""
    wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    ws = wb[sheet] if sheet in wb.sheetnames else wb[wb.sheetnames[0]]
    mapping = {}
    dups = 0
    header = None
    e_i = n_i = None
    for row in ws.iter_rows(values_only=True):
        if header is None:
            header = [str(c).strip() if c is not None else "" for c in row]
            try:
                e_i = header.index(entry_col)
                n_i = header.index(name_col)
            except ValueError:
                raise SystemExit(
                    f"UniProt xlsx missing required column(s). "
                    f"Wanted '{entry_col}' and '{name_col}'; got {header}")
            continue
        if row is None or len(row) <= max(e_i, n_i):
            continue
        entry = str(row[e_i]).strip() if row[e_i] is not None else ""
        if not entry:
            continue
        name = str(row[n_i]).strip() if row[n_i] is not None else ""
        if entry in mapping:
            dups += 1
            if not mapping[entry] and name:
                mapping[entry] = name
        else:
            mapping[entry] = name
    wb.close()
    return mapping, dups


def resolve_uniprot(accession: str, mapping: dict):
    """(label, status). Exact first, then canonical fallback for isoforms."""
    acc = str(accession or "").strip()
    if not acc:
        return "", "missing"
    if acc in mapping:
        return mapping[acc], "exact"
    m = ISOFORM_SUFFIX_RE.fullmatch(acc)
    if m:
        canonical = m.group("canonical")
        if canonical in mapping:
            return mapping[canonical], "canonical_fallback"
    return "", "missing"


# ---------------------------------------------------------------------------
# besthits.tsv reader
# ---------------------------------------------------------------------------
def read_besthits(path: Path):
    """Return {(part, db): {qseqid: row}} plus the set of accessions per db."""
    by_part_db = defaultdict(dict)
    accs_by_db = defaultdict(set)
    with path.open("r", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (row["part"], row["db"])
            by_part_db[key][row["qseqid"]] = row
            if row.get("accession"):
                accs_by_db[row["db"]].add(row["accession"])
    return by_part_db, accs_by_db


# ---------------------------------------------------------------------------
# xlsx annotation
# ---------------------------------------------------------------------------
def _find_col(headers, name):
    for i, h in enumerate(headers, start=1):
        if h is not None and str(h).strip() == name:
            return i
    return None


def annotate_workbook(xlsx_in: Path, xlsx_out: Path, key_col_name: str,
                      labels_by_key: dict):
    """Write one column per DB into every sheet that carries the key column."""
    wb = openpyxl.load_workbook(xlsx_in)
    rows_total = 0
    filled = Counter()
    for sh in wb.sheetnames:
        ws = wb[sh]
        headers = [c.value for c in ws[1]]
        key_col = _find_col(headers, key_col_name)
        if key_col is None:
            continue
        out_cols = {}
        next_col = len(headers) + 1
        for db in DB_TAGS:
            c = _find_col(headers, db)
            if c is None:
                c = next_col
                ws.cell(row=1, column=c, value=db)
                next_col += 1
            out_cols[db] = c
        for r in range(2, ws.max_row + 1):
            key = ws.cell(row=r, column=key_col).value
            if key is None:
                continue
            key = str(key).strip()
            if not key:
                continue
            rows_total += 1
            per_db = labels_by_key.get(key, {})
            for db in DB_TAGS:
                val = per_db.get(db, "")
                if val:
                    filled[db] += 1
                ws.cell(row=r, column=out_cols[db], value=val)
    xlsx_out.parent.mkdir(parents=True, exist_ok=True)
    wb.save(xlsx_out)
    return rows_total, filled


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--besthits",           required=True, type=Path)
    ap.add_argument("--pos-xlsx",           required=True, type=Path)
    ap.add_argument("--neg-xlsx",           required=True, type=Path)
    ap.add_argument("--conoserver-fasta",   required=True, type=Path)
    ap.add_argument("--conoserver-field",   type=int, default=6)
    ap.add_argument("--ncbi-feature-table", required=True, type=Path)
    ap.add_argument("--uniprot-xlsx",       required=True, type=Path)
    ap.add_argument("--uniprot-sheet",      default="Sheet0")
    ap.add_argument("--uniprot-entry-col",  default="Entry")
    ap.add_argument("--uniprot-name-col",   default="Protein names")
    ap.add_argument("--key-col",            default="seq_id")
    ap.add_argument("--out-pos-xlsx",       required=True, type=Path)
    ap.add_argument("--out-neg-xlsx",       required=True, type=Path)
    ap.add_argument("--out-labels",         required=True, type=Path)
    ap.add_argument("--out-audit",          required=True, type=Path)
    args = ap.parse_args()

    by_part_db, accs_by_db = read_besthits(args.besthits)

    # ---- Load the three label sources -----------------------------------
    cs_idx, cs_dups = load_conoserver_labels(args.conoserver_fasta,
                                             args.conoserver_field)
    ncbi_idx = load_ncbi_labels(args.ncbi_feature_table, accs_by_db["ncbi"])
    up_idx, up_dups = load_uniprot_labels(args.uniprot_xlsx, args.uniprot_sheet,
                                          args.uniprot_entry_col,
                                          args.uniprot_name_col)

    # ---- Resolve labels per (part, qseqid) -------------------------------
    labels = defaultdict(lambda: defaultdict(dict))   # part -> key -> db -> label
    status_counts = Counter()
    missing_accs = Counter()

    for (part, db), rows in by_part_db.items():
        for q, row in rows.items():
            if db == "conoserver":
                # Collapse the distinct superfamily labels of ALL hits.
                vals, seen = [], set()
                for sseqid in (row.get("all_sseqids") or "").split(";"):
                    acc = sseqid.split("|", 1)[0].strip()
                    if not acc:
                        continue
                    if acc not in cs_idx:
                        missing_accs[f"conoserver:{acc}"] += 1
                        continue
                    v = (cs_idx.get(acc) or "").strip()
                    if v and v not in seen:
                        vals.append(v)
                        seen.add(v)
                label = ", ".join(vals)
                status_counts[f"conoserver:{'ok' if label else 'empty'}"] += 1
            elif db == "ncbi":
                acc = row.get("accession", "")
                if acc and acc not in ncbi_idx:
                    missing_accs[f"ncbi:{acc}"] += 1
                label = (ncbi_idx.get(acc) or "").strip()
                status_counts[f"ncbi:{'ok' if label else 'empty'}"] += 1
            else:  # uniprot
                label, st = resolve_uniprot(row.get("accession", ""), up_idx)
                status_counts[f"uniprot:{st}"] += 1
                if st == "missing" and row.get("accession"):
                    missing_accs[f"uniprot:{row['accession']}"] += 1
            if label:
                labels[part][q][db] = label

    # ---- Write annotated workbooks --------------------------------------
    pos_total, pos_filled = annotate_workbook(
        args.pos_xlsx, args.out_pos_xlsx, args.key_col, labels["positivo"])
    neg_total, neg_filled = annotate_workbook(
        args.neg_xlsx, args.out_neg_xlsx, args.key_col, labels["negativo"])

    # ---- Tidy labels.tsv -------------------------------------------------
    args.out_labels.parent.mkdir(parents=True, exist_ok=True)
    with args.out_labels.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["seq_id", "part"] + list(DB_TAGS))
        for part in ("positivo", "negativo"):
            for q in sorted(labels[part]):
                per_db = labels[part][q]
                w.writerow([q, part] + [per_db.get(db, "") for db in DB_TAGS])

    # ---- Audit -----------------------------------------------------------
    args.out_audit.parent.mkdir(parents=True, exist_ok=True)
    with args.out_audit.open("w", encoding="utf-8", newline="") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"conoserver_db_entries\t{len(cs_idx)}\n")
        fh.write(f"conoserver_dup_accessions\t{len(cs_dups)}\n")
        fh.write(f"ncbi_blocks_resolved\t{len(ncbi_idx)}\n")
        fh.write(f"ncbi_accessions_queried\t{len(accs_by_db['ncbi'])}\n")
        fh.write(f"uniprot_db_entries\t{len(up_idx)}\n")
        fh.write(f"uniprot_dup_entries\t{up_dups}\n")
        fh.write(f"rows_positivo\t{pos_total}\n")
        fh.write(f"rows_negativo\t{neg_total}\n")
        for db in DB_TAGS:
            fh.write(f"filled_positivo_{db}\t{pos_filled[db]}\n")
        for db in DB_TAGS:
            fh.write(f"filled_negativo_{db}\t{neg_filled[db]}\n")
        for k in sorted(status_counts):
            fh.write(f"status_{k.replace(':', '_')}\t{status_counts[k]}\n")
        fh.write(f"distinct_missing_accessions\t{len(missing_accs)}\n")

    print(f"[annotate_all_labels] pos_rows={pos_total} neg_rows={neg_total} "
          f"filled_pos={dict(pos_filled)} filled_neg={dict(neg_filled)}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
