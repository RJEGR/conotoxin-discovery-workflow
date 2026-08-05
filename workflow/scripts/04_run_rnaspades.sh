#!/usr/bin/env bash
# =============================================================================
# 04_run_rnaspades.sh  -  rnaSPAdes de novo assembly (SE/PE)
# =============================================================================
# Runs natively from a conda env (no container): SPAdes ships native builds
# for linux-64 and osx-64, and on Apple Silicon it runs under Rosetta without
# intervention.
#
# CANONICAL CLI
#   PE:  rnaspades.py -1 R1 -2 R2 -o OUTDIR -t THREADS -m MEM_GB
#   SE:  rnaspades.py -s R1      -o OUTDIR -t THREADS -m MEM_GB
#
# rnaSPAdes emits OUTDIR/transcripts.fasta (plus hard_filtered_transcripts.fasta
# from 3.15 onward). We take `transcripts.fasta`: filtering is Stage 01's job
# via TransRate, and applying SPAdes' own filter first would compose two
# different contig-selection criteria — one of them assembler-specific — which
# would make the three assemblers no longer comparable.
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

OUTDIR="$(dirname "$CONTIGS_OUT")/rnaspades_work"
mkdir -p "$OUTDIR"

echo "[rnaspades] sample=$SAMPLE_ID layout=$LAYOUT threads=$THREADS mem=${MEM_GB}G"
VERSION_STR="$(rnaspades.py --version 2>&1 | head -1 | tr -d '\n')"
echo "[rnaspades] version=$VERSION_STR"

CMD=(rnaspades.py -o "$OUTDIR" -t "$THREADS" -m "$MEM_GB")
if [ "$LAYOUT" = "PE" ]; then
    [ -n "$R2" ] || { echo "[rnaspades] ERROR: PE but no R2" >&2; exit 2; }
    CMD+=(-1 "$R1" -2 "$R2")
else
    CMD+=(-s "$R1")
fi

echo "[rnaspades] cmd: ${CMD[*]}"
t_start=$(date +%s)
"${CMD[@]}"
t_end=$(date +%s)
ELAPSED=$(( t_end - t_start ))

NATIVE_OUT="${OUTDIR}/transcripts.fasta"
if [ ! -s "$NATIVE_OUT" ]; then
    echo "[rnaspades] ERROR: no transcripts.fasta in $OUTDIR" >&2
    exit 3
fi
mv "$NATIVE_OUT" "$CONTIGS_OUT"

# SPAdes leaves a large K-mer scaffold tree behind; keep only the log for
# provenance. On a 40-sample run the intermediates are tens of GB.
if [ -f "${OUTDIR}/spades.log" ]; then
    mv "${OUTDIR}/spades.log" "$(dirname "$CONTIGS_OUT")/rnaspades.log"
fi
rm -rf "$OUTDIR"

read -r N_CONTIGS N_BASES <<< "$(fasta_stats "$CONTIGS_OUT")"
{
    printf 'metric\tvalue\n'
    printf 'assembler\trnaspades\n'
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

echo "[rnaspades] done: contigs=$N_CONTIGS bases=$N_BASES elapsed=${ELAPSED}s"
