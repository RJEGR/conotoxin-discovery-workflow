#!/usr/bin/env bash
# =============================================================================
# 06_run_transrate.sh  -  reference-free contig filtering, and the
#                          filtered_contigs.fasta contract
# =============================================================================
# Produces `filtered_contigs.fasta`, the SINGLE input Stage 02 consumes. What
# fills it depends on the library layout:
#
#   PAIRED-END -> TransRate read-metrics mode. TransRate maps the reads back
#     onto the assembly, scores every contig on the resulting evidence, and
#     writes `good.<stem>.fa` containing the well-supported subset. That file
#     becomes filtered_contigs.fasta and `filter_method` is recorded as
#     `transrate`.
#
#   SINGLE-END -> no filtering. TransRate's read-metrics mode requires BOTH
#     --left and --right; there is no documented single-end path, and the
#     reference mode needs a reference transcriptome that does not exist for
#     real data. The raw assembly is copied through unchanged and
#     `filter_method` is recorded as `none`.
#
#     CONSEQUENCE (decision 2026-07-30): SE and PE contig sets are NOT
#     directly comparable — the SE sets retain low-support contigs that the PE
#     sets have removed. Any downstream comparison must stratify on
#     `filter_method`, which is why it is a manifest column rather than a
#     footnote.
#
# TransRate 1.0.3 predates transparent gzip input, so reads are decompressed
# into a scratch directory first.
#
# USAGE (positional)
#   $1 layout(SE|PE) $2 contigs $3 r1 $4 r2("" if SE)
#   $5 filtered_out $6 summary_out $7 transrate_dir
#   $8 threads $9 sample_id $10 assembler $11 enabled(true|false)
#   $12 sanitize(true|false)
#   $13 runtime $14 image $15 platform $16 cache_dir $17 extra_binds
# =============================================================================

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib_container_paths.sh
source "${HERE}/lib_container_paths.sh"

LAYOUT="${1:?}"; CONTIGS="${2:?}"; R1="${3:?}"; R2="${4:-}"
FILTERED_OUT="${5:?}"; SUMMARY="${6:?}"; TR_DIR="${7:?}"
THREADS="${8:?}"; SAMPLE_ID="${9:?}"; ASSEMBLER="${10:?}"
ENABLED="${11:-true}"; SANITIZE="${12:-true}"
RUNTIME="${13:?}"; IMAGE="${14:?}"; PLATFORM="${15:-}"
CACHE_DIR="${16:-resources/containers}"; EXTRA_BINDS="${17:-}"

mkdir -p "$TR_DIR" "$(dirname "$FILTERED_OUT")"

read -r N_IN _ <<< "$(fasta_stats "$CONTIGS")"
echo "[transrate] sample=$SAMPLE_ID asm=$ASSEMBLER layout=$LAYOUT contigs_in=$N_IN"

emit_summary() {
    local method="$1" score="$2" elapsed="$3" note="$4"
    local n_out n_bases
    read -r n_out n_bases <<< "$(fasta_stats "$FILTERED_OUT")"
    local pct="0"
    [ "$N_IN" -gt 0 ] && pct=$(awk -v a="$n_out" -v b="$N_IN" \
                                'BEGIN{printf "%.2f", (a/b)*100}')
    {
        printf 'metric\tvalue\n'
        printf 'sample_id\t%s\n'          "$SAMPLE_ID"
        printf 'assembler\t%s\n'          "$ASSEMBLER"
        printf 'layout\t%s\n'             "$LAYOUT"
        printf 'filter_method\t%s\n'      "$method"
        printf 'transrate_score\t%s\n'    "$score"
        printf 'n_contigs_in\t%s\n'       "$N_IN"
        printf 'n_contigs_filtered\t%s\n' "$n_out"
        printf 'n_bases_filtered\t%s\n'   "$n_bases"
        printf 'pct_contigs_retained\t%s\n' "$pct"
        printf 'filtered_fasta\t%s\n'     "$FILTERED_OUT"
        printf 'elapsed_seconds\t%s\n'    "$elapsed"
        printf 'notes\t%s\n'              "$note"
    } > "$SUMMARY"
    echo "[transrate] method=$method contigs $N_IN -> $n_out (${pct}%)"
}

# ---- Cases that bypass TransRate ------------------------------------------
if [ "$ENABLED" != "true" ]; then
    cp "$CONTIGS" "$FILTERED_OUT"
    emit_summary "none" "NA" "0" "transrate_disabled_in_config"
    exit 0
fi
if [ "$LAYOUT" = "SE" ]; then
    cp "$CONTIGS" "$FILTERED_OUT"
    emit_summary "none" "NA" "0" "single_end_transrate_requires_paired_reads"
    exit 0
fi
if [ "$N_IN" -eq 0 ]; then
    cp "$CONTIGS" "$FILTERED_OUT"
    emit_summary "none" "NA" "0" "empty_assembly"
    exit 0
fi

# ---- Paired-end: run TransRate --------------------------------------------
# DISK BUDGET. TransRate 1.0.3 cannot read gzipped FASTQ, so the reads must be
# materialised uncompressed. That is expensive: ~350 bytes per record means a
# 9.3 M-pair library needs ~7 GB for both mates.
#
# An earlier version of this script staged the reads FOUR times (decompress to
# scratch -> seqkit sana -> seqkit pair -> copy into the container-visible
# directory), i.e. ~26 GB of temporaries, and ran out of space on a 20 GB
# volume. Everything below writes each mate exactly ONCE, straight into the
# directory the container reads from.
STEM="$(basename "$TR_DIR")"
STAGE_DIR="${TR_DIR}/_input"
mkdir -p "$STAGE_DIR"

# Always clean the staged reads, including on failure — leaving 7 GB behind
# after a crash is how the next run runs out of space too.
cleanup() { rm -rf "$STAGE_DIR"; }
trap cleanup EXIT

# ---- Preflight: is there room? --------------------------------------------
# Estimate the uncompressed size as 4x the gzipped size (typical for FASTQ)
# and refuse to start if the filesystem cannot hold both mates plus 20% slack.
# Failing here with a clear number beats failing halfway through a pigz write.
NEED_KB=$(( ( $(du -k "$R1" | cut -f1) + $(du -k "$R2" | cut -f1) ) * 4 * 12 / 10 ))
AVAIL_KB=$(df -k "$STAGE_DIR" | awk 'NR==2{print $4}')
if [ "$AVAIL_KB" -lt "$NEED_KB" ]; then
    echo "[transrate] ERROR: need ~$(( NEED_KB / 1024 )) MB to decompress the" \
         "reads but only $(( AVAIL_KB / 1024 )) MB are free on the volume" \
         "holding $STAGE_DIR." >&2
    echo "[transrate] Free space, or set assembly.transrate.enabled=false to" \
         "skip filtering (filtered_contigs.fasta then carries the raw" \
         "assembly and filter_method is recorded as 'none')." >&2
    cp "$CONTIGS" "$FILTERED_OUT"
    emit_summary "none" "NA" "0" "insufficient_disk_for_transrate"
    exit 0
fi

# ---- Materialise the inputs (one write each) ------------------------------
cp "$CONTIGS" "${STAGE_DIR}/${STEM}.fa"

decompress_to() {
    local src="$1" dst="$2"
    case "$src" in
        *.gz) if command -v pigz >/dev/null 2>&1; then
                  pigz -dc "$src" > "$dst"
              else
                  gzip -dc "$src" > "$dst"
              fi ;;
        *)    cp "$src" "$dst" ;;
    esac
}
decompress_to "$R1" "${STAGE_DIR}/R1.fq"
decompress_to "$R2" "${STAGE_DIR}/R2.fq"

# Read-pair sanitisation. Off by default for this pipeline: Stage 00 runs fastp
# with --in1/--in2, which emits synchronised mates by construction, so there is
# nothing to repair. (The benchmark pipeline needs it because ART occasionally
# simulates out-of-sync mates.) When enabled, the repair writes to a temporary
# and replaces the original in place, so the peak footprint grows by one mate
# rather than by three.
if [ "$SANITIZE" = "true" ] && command -v seqkit >/dev/null 2>&1; then
    echo "[transrate] sanitising read pairing with seqkit (in place)"
    for mate in R1 R2; do
        if seqkit sana "${STAGE_DIR}/${mate}.fq" \
                       -o "${STAGE_DIR}/${mate}.sana.fq" 2>/dev/null; then
            mv "${STAGE_DIR}/${mate}.sana.fq" "${STAGE_DIR}/${mate}.fq"
        else
            rm -f "${STAGE_DIR}/${mate}.sana.fq"
        fi
    done
fi

QUERY_IN="$(to_container_path "${STAGE_DIR}/${STEM}.fa" "$EXTRA_BINDS")"
R1_IN="$(to_container_path "${STAGE_DIR}/R1.fq" "$EXTRA_BINDS")"
R2_IN="$(to_container_path "${STAGE_DIR}/R2.fq" "$EXTRA_BINDS")"
OUT_IN="$(to_container_path "$(dirname "$TR_DIR")" "$EXTRA_BINDS")"

TR_CMD=(transrate --assembly "$QUERY_IN"
        --left "$R1_IN" --right "$R2_IN"
        --output "$OUT_IN" --threads "$THREADS")

echo "[transrate] cmd: ${TR_CMD[*]}"
t_start=$(date +%s)
set +e
bash "${HERE}/container_run.sh" "$RUNTIME" "$IMAGE" "$PLATFORM" \
     "$CACHE_DIR" "$EXTRA_BINDS" -- "${TR_CMD[@]}"
RC=$?
set -e
t_end=$(date +%s)
ELAPSED=$(( t_end - t_start ))
echo "[transrate] container exit=$RC elapsed=${ELAPSED}s"

# ---- Environment failures must NOT degrade silently ------------------------
# container_run.sh reserves exit 10 (no usable runtime) and 11 (Docker daemon
# unreachable) for environment problems, and 12 for a failed image pull. Those
# are setup faults the operator fixes in seconds — not statements about the
# data — so they abort the rule.
#
# Treating them like a TransRate failure (i.e. passing the raw contigs through)
# was the earlier behaviour and it actively hid the problem: on 2026-08-03 a
# stopped Docker Desktop produced two "successful" assemblies whose contigs
# were entirely unfiltered, recorded only as a note nobody reads until the
# annotation looks odd. Filtered and unfiltered contig sets are not
# comparable, so producing them by accident is worse than stopping.
#
# Every other exit code still degrades gracefully below: those mean TransRate
# ran but could not score THIS assembly, which is a legitimate data outcome
# and should not abort a multi-sample batch.
case "$RC" in
    10|11|12)
        echo "[transrate] ERROR: container runtime unavailable (exit=$RC)." >&2
        echo "[transrate] This is an environment fault, not a data problem, so" \
             "the rule fails instead of silently emitting unfiltered contigs." >&2
        echo "[transrate] macOS: start Docker Desktop and check" \
             "\`docker info\`. HPC: load the apptainer/singularity module." >&2
        echo "[transrate] To proceed deliberately without filtering, set" \
             "assembly.transrate.enabled=false in config.yaml." >&2
        exit "$RC"
        ;;
esac

GOOD_FA="${TR_DIR}/good.${STEM}.fa"
[ -s "$GOOD_FA" ] || GOOD_FA=$( { find "$TR_DIR" -maxdepth 2 -type f \
                                   -name 'good.*.fa' -print -quit; } \
                                2>/dev/null || true )

SCORE="NA"
# TransRate writes assemblies.csv at the ROOT of --output (i.e. the assembler
# directory), while contigs.csv and good.*.fa go into <output>/<stem>/. We pass
# --output "$(dirname "$TR_DIR")", so the summary CSV is one level ABOVE
# $TR_DIR. Looking for it inside $TR_DIR silently yielded score=NA on every
# run until 2026-08-04.
ASSEMBLIES_CSV="$(dirname "$TR_DIR")/assemblies.csv"
[ -s "$ASSEMBLIES_CSV" ] || ASSEMBLIES_CSV="${TR_DIR}/assemblies.csv"
if [ -s "$ASSEMBLIES_CSV" ]; then
    SCORE=$(python - "$ASSEMBLIES_CSV" << 'PY'
import csv, sys
try:
    with open(sys.argv[1], newline="") as fh:
        row = next(csv.DictReader(fh), {})
    print(row.get("score", "NA") or "NA")
except Exception:
    print("NA")
PY
)
fi

rm -rf "$STAGE_DIR"
trap - EXIT

if [ -s "${GOOD_FA:-}" ]; then
    cp "$GOOD_FA" "$FILTERED_OUT"
    emit_summary "transrate" "$SCORE" "$ELAPSED" ""
else
    # TransRate ran but could not score this assembly (too few reads mapped,
    # a degenerate contig set, an internal SNAP error...). That is a data
    # outcome, not an environment fault, so the sample stays in the run rather
    # than aborting a multi-sample batch. The note and `filter_method = none`
    # make the degradation explicit in the manifest — downstream comparisons
    # must stratify on it, because unfiltered contig sets are not comparable
    # with filtered ones.
    echo "[transrate] WARN: no good.*.fa produced (exit=$RC); passing raw contigs through" >&2
    cp "$CONTIGS" "$FILTERED_OUT"
    emit_summary "none" "$SCORE" "$ELAPSED" "transrate_failed_exit_${RC}"
fi
