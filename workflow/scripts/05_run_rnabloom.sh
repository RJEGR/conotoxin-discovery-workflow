#!/usr/bin/env bash
# =============================================================================
# 05_run_rnabloom.sh  -  RNA-Bloom de novo assembly (SE/PE)
# =============================================================================
# Runs natively from a conda env: RNA-Bloom is pure Java, so Bioconda ships a
# noarch build that works on osx-arm64, osx-64 and linux-64 alike.
#
# CANONICAL CLI
#   PE:  rnabloom -left R1 -right R2 -revcomp-right -t N -mem G -outdir DIR
#   SE:  rnabloom -sef R1                           -t N -mem G -outdir DIR
#
# `-revcomp-right` is REQUIRED for standard Illumina paired-end data. Illumina
# PE reads are in FR orientation (R2 sequenced from the opposite strand), and
# RNA-Bloom expects both mates on the same strand. Omitting the flag does not
# error — it silently fragments the assembly, which is the worst kind of bug
# because it looks like a biological result.
#
# USAGE (positional)
#   $1 layout(SE|PE) $2 r1 $3 r2("" if SE) $4 contigs_out $5 summary_out
#   $6 mem_gb $7 threads $8 sample_id
# =============================================================================

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib_container_paths.sh
source "${HERE}/lib_container_paths.sh"

LAYOUT="${1:?}"; R1="${2:?}"; R2="${3:-}"
CONTIGS_OUT="${4:?}"; SUMMARY="${5:?}"
MEM_GB="${6:?}"; THREADS="${7:?}"; SAMPLE_ID="${8:?}"

OUTDIR="$(dirname "$CONTIGS_OUT")/rnabloom_work"
mkdir -p "$OUTDIR"

echo "[rnabloom] sample=$SAMPLE_ID layout=$LAYOUT threads=$THREADS mem=${MEM_GB}G"
VERSION_STR="$(rnabloom --version 2>&1 | head -1 | tr -d '\n' || echo unknown)"
echo "[rnabloom] version=$VERSION_STR"

CMD=(rnabloom -t "$THREADS" -mem "$MEM_GB" -outdir "$OUTDIR")
if [ "$LAYOUT" = "PE" ]; then
    [ -n "$R2" ] || { echo "[rnabloom] ERROR: PE but no R2" >&2; exit 2; }
    CMD+=(-left "$R1" -right "$R2" -revcomp-right)
else
    CMD+=(-sef "$R1")
fi

echo "[rnabloom] cmd: ${CMD[*]}"
t_start=$(date +%s)
"${CMD[@]}"
t_end=$(date +%s)
ELAPSED=$(( t_end - t_start ))

NATIVE_OUT="${OUTDIR}/rnabloom.transcripts.fa"
if [ ! -s "$NATIVE_OUT" ]; then
    NATIVE_OUT=$( { find "$OUTDIR" -maxdepth 2 -type f \
                     \( -name 'rnabloom.transcripts.fa' \
                        -o -name '*transcripts*.fa' \
                        -o -name '*transcripts*.fasta' \) -print -quit; } \
                  2>/dev/null || true )
fi
if [ ! -s "${NATIVE_OUT:-}" ]; then
    echo "[rnabloom] ERROR: no transcripts fasta under $OUTDIR" >&2
    exit 3
fi
mv "$NATIVE_OUT" "$CONTIGS_OUT"

# RNA-Bloom leaves Bloom-filter binaries behind; they are large and worthless
# once the assembly exists.
rm -rf "$OUTDIR"

read -r N_CONTIGS N_BASES <<< "$(fasta_stats "$CONTIGS_OUT")"
{
    printf 'metric\tvalue\n'
    printf 'assembler\trnabloom\n'
    printf 'assembler_version\t%s\n' "$VERSION_STR"
    printf 'container_image\t%s\n'   "native-conda"
    printf 'sample_id\t%s\n'         "$SAMPLE_ID"
    printf 'layout\t%s\n'            "$LAYOUT"
    printf 'threads\t%s\n'           "$THREADS"
    printf 'max_memory_gb\t%s\n'     "$MEM_GB"
    printf 'elapsed_seconds\t%s\n'   "$ELAPSED"
    printf 'n_contigs\t%s\n'         "$N_CONTIGS"
    printf 'n_bases\t%s\n'           "$N_BASES"
    printf 'contigs_fasta\t%s\n'     "$CONTIGS_OUT"
} > "$SUMMARY"

echo "[rnabloom] done: contigs=$N_CONTIGS bases=$N_BASES elapsed=${ELAPSED}s"
