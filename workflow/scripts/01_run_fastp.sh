#!/usr/bin/env bash
# =============================================================================
# 01_run_fastp.sh  -  Stage 00 wrapper for fastp (single- and paired-end)
# =============================================================================
# One wrapper handles both layouts: the caller passes an empty string for
# `in_r2` to request single-end mode. Keeping the branch here rather than in
# two Snakemake rules means the emitted summary schema, the log format and the
# post-run verification are identical for SE and PE — the manifest builder and
# every downstream consumer then need no layout-specific code.
#
# USAGE (positional)
#    $1  layout            SE | PE
#    $2  in_r1             input FASTQ (R1, or the only file when SE)
#    $3  in_r2             input FASTQ R2, or "" when SE
#    $4  out_r1            cleaned R1
#    $5  out_r2            cleaned R2, or "" when SE
#    $6  html              QC report (HTML)
#    $7  json              QC report (JSON)
#    $8  summary_tsv       per-sample metrics for the stage manifest
#    $9  threads
#   $10  sample_id
#   $11..  extra fastp flags, already expanded by the rule
#
# OUTPUT `summary_tsv` schema (key/value, one pair per line):
#   sample_id, layout, fastp_version, threads,
#   reads_before, reads_after, bases_before, bases_after,
#   q20_rate_after, q30_rate_after, gc_content_after, duplication_rate,
#   read1_mean_length_after, read2_mean_length_after,
#   pct_reads_passed, elapsed_seconds
#
# Metrics come from fastp's own JSON (authoritative, computed in-stream). The
# emitted FASTQs are additionally re-counted with seqkit to detect truncated
# writes: fastp can die mid-write and still leave a syntactically valid but
# incomplete gzip member, which would silently shrink an assembly.
# =============================================================================

set -euo pipefail

LAYOUT="${1:?layout required}"
IN_R1="${2:?in_r1 required}"
IN_R2="${3:-}"
OUT_R1="${4:?out_r1 required}"
OUT_R2="${5:-}"
HTML="${6:?html required}"
JSON="${7:?json required}"
SUMMARY="${8:?summary required}"
THREADS="${9:-4}"
SAMPLE_ID="${10:?sample_id required}"
shift 10
# `set -u` makes `"${EXTRA[@]}"` an error on an unset/empty array in older
# bash, so guard every expansion with the `+` form.
EXTRA=()
[ "$#" -gt 0 ] && EXTRA=("$@")

mkdir -p "$(dirname "$OUT_R1")" "$(dirname "$HTML")" "$(dirname "$SUMMARY")"

echo "[fastp] sample=$SAMPLE_ID layout=$LAYOUT threads=$THREADS"
echo "[fastp] in_r1=$IN_R1"
[ -n "$IN_R2" ] && echo "[fastp] in_r2=$IN_R2"

FASTP_VERSION="$(fastp --version 2>&1 | head -1 | tr -d '\n')"
echo "[fastp] version=$FASTP_VERSION"

# ---- Build the command ----------------------------------------------------
CMD=(fastp
        --in1 "$IN_R1"
        --out1 "$OUT_R1"
        --html "$HTML"
        --json "$JSON"
        --thread "$THREADS"
        --report_title "fastp report - ${SAMPLE_ID}")

if [ "$LAYOUT" = "PE" ]; then
    if [ -z "$IN_R2" ] || [ -z "$OUT_R2" ]; then
        echo "[fastp] ERROR: layout=PE but in_r2/out_r2 not supplied" >&2
        exit 2
    fi
    CMD+=(--in2 "$IN_R2" --out2 "$OUT_R2")
fi

CMD+=(${EXTRA[@]+"${EXTRA[@]}"})

echo "[fastp] cmd: ${CMD[*]}"

t_start=$(date +%s)
"${CMD[@]}"
t_end=$(date +%s)
ELAPSED=$(( t_end - t_start ))

# ---- Verify the emitted files ---------------------------------------------
# A zero-byte or unreadable output means fastp reported success but wrote
# nothing usable; fail loudly here rather than letting an empty assembly
# propagate three stages downstream.
verify_fastq() {
    local f="$1" label="$2"
    if [ ! -s "$f" ]; then
        echo "[fastp] ERROR: $label is missing or empty: $f" >&2
        exit 3
    fi
    if ! gzip -t "$f" 2>/dev/null; then
        echo "[fastp] ERROR: $label is not a valid gzip stream (truncated write?): $f" >&2
        exit 3
    fi
}
verify_fastq "$OUT_R1" "out_r1"
[ "$LAYOUT" = "PE" ] && verify_fastq "$OUT_R2" "out_r2"

# Independent read count from the emitted file (fastp's JSON says what it
# *intended* to write; seqkit says what is actually on disk).
SEQKIT_R1=$(seqkit stats -T "$OUT_R1" 2>/dev/null | awk 'NR==2{print $4}')
SEQKIT_R1="${SEQKIT_R1:-NA}"

# ---- Extract metrics from fastp's JSON ------------------------------------
# Parsed with python rather than jq: python is already guaranteed by the env
# stack, jq is not, and this keeps the env one package smaller.
python - "$JSON" "$SUMMARY" "$SAMPLE_ID" "$LAYOUT" "$FASTP_VERSION" \
        "$THREADS" "$ELAPSED" "$SEQKIT_R1" "$OUT_R1" "$OUT_R2" << 'PY'
import json, sys

(json_path, out_path, sample_id, layout, version,
 threads, elapsed, seqkit_r1, out_r1, out_r2) = sys.argv[1:11]

with open(json_path) as fh:
    d = json.load(fh)

before = d.get("summary", {}).get("before_filtering", {})
after  = d.get("summary", {}).get("after_filtering", {})
filt   = d.get("filtering_result", {})

reads_before = int(before.get("total_reads", 0) or 0)
reads_after  = int(after.get("total_reads", 0) or 0)
pct = round(reads_after / reads_before * 100, 2) if reads_before else 0.0

rows = [
    ("sample_id",               sample_id),
    ("layout",                  layout),
    ("fastp_version",           version),
    ("threads",                 threads),
    ("reads_before",            reads_before),
    ("reads_after",             reads_after),
    ("bases_before",            before.get("total_bases", 0)),
    ("bases_after",             after.get("total_bases", 0)),
    ("q20_rate_after",          after.get("q20_rate", "")),
    ("q30_rate_after",          after.get("q30_rate", "")),
    ("gc_content_after",        after.get("gc_content", "")),
    ("read1_mean_length_after", after.get("read1_mean_length", "")),
    ("read2_mean_length_after", after.get("read2_mean_length", "") if layout == "PE" else ""),
    ("duplication_rate",        d.get("duplication", {}).get("rate", "")),
    ("reads_passed_filter",     filt.get("passed_filter_reads", "")),
    ("reads_low_quality",       filt.get("low_quality_reads", "")),
    ("reads_too_many_N",        filt.get("too_many_N_reads", "")),
    ("reads_too_short",         filt.get("too_short_reads", "")),
    ("pct_reads_passed",        pct),
    ("seqkit_reads_out_r1",     seqkit_r1),
    ("clean_r1",                out_r1),
    ("clean_r2",                out_r2),
    ("elapsed_seconds",         elapsed),
]

with open(out_path, "w") as fh:
    fh.write("metric\tvalue\n")
    for k, v in rows:
        fh.write(f"{k}\t{v}\n")

print(f"[fastp] reads {reads_before:,} -> {reads_after:,} ({pct}% retained)")
PY

echo "[fastp] done in ${ELAPSED}s -> $OUT_R1"
