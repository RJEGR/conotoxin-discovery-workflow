#!/usr/bin/env bash
# =============================================================================
# 03_run_trinity.sh  -  Trinity de novo assembly (SE/PE) in a container
# =============================================================================
# Trinity runs inside trinityrnaseq/trinityrnaseq:2.15.x, invoked through
# container_run.sh so the same code path works with Docker on macOS and
# Apptainer on SLURM.
#
# WHY 2.15 AND NOT THE BENCHMARK'S 2.6.6. The simulation benchmark pins
# biocontainers 2.6.6 (2018) because Bioconda's Trinity was unsolvable on
# Apple Silicon at the time. For real, deep libraries the modern release
# assembles better, handles memory more predictably, and supports --single
# cleanly. Assemblies from this pipeline are therefore NOT byte-comparable
# with the benchmark's — do not merge the two in one table.
#
# CANONICAL CLI
#   PE:  Trinity --seqType fq --left R1 --right R2 --max_memory NG --CPU N \
#                --output <dir containing 'trinity'> --full_cleanup
#   SE:  Trinity --seqType fq --single R1  ...
#
# USAGE (positional)
#   $1 layout(SE|PE) $2 r1 $3 r2("" if SE) $4 contigs_out $5 summary_out
#   $6 mem_gb $7 threads $8 sample_id $9 no_normalize(true|false)
#   $10 runtime $11 image $12 platform $13 cache_dir $14 extra_binds
# =============================================================================

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib_container_paths.sh
source "${HERE}/lib_container_paths.sh"

LAYOUT="${1:?}"; R1="${2:?}"; R2="${3:-}"
CONTIGS_OUT="${4:?}"; SUMMARY="${5:?}"
MEM_GB="${6:?}"; THREADS="${7:?}"; SAMPLE_ID="${8:?}"
NO_NORMALIZE="${9:-false}"
RUNTIME="${10:?}"; IMAGE="${11:?}"; PLATFORM="${12:-}"
CACHE_DIR="${13:-resources/containers}"; EXTRA_BINDS="${14:-}"

OUTDIR_HOST="$(dirname "$CONTIGS_OUT")"
mkdir -p "$OUTDIR_HOST"

echo "[trinity] sample=$SAMPLE_ID layout=$LAYOUT threads=$THREADS mem=${MEM_GB}G"
echo "[trinity] image=$IMAGE runtime=$RUNTIME"

R1_IN="$(to_container_path "$R1" "$EXTRA_BINDS")"
OUTDIR_IN="$(to_container_path "$OUTDIR_HOST" "$EXTRA_BINDS")"

# Trinity refuses any --output path that does not contain the substring
# 'trinity' (case-insensitive), so the leaf is fixed. With --full_cleanup the
# assembly lands next to that directory as <leaf>.Trinity.fasta.
TRINITY_OUT_IN="${OUTDIR_IN}/trinity_out_dir"
NATIVE_OUT_HOST="${OUTDIR_HOST}/trinity_out_dir.Trinity.fasta"

TRINITY_CMD=(Trinity --seqType fq --max_memory "${MEM_GB}G" --CPU "$THREADS"
             --output "$TRINITY_OUT_IN" --full_cleanup)

if [ "$LAYOUT" = "PE" ]; then
    [ -n "$R2" ] || { echo "[trinity] ERROR: PE but no R2" >&2; exit 2; }
    R2_IN="$(to_container_path "$R2" "$EXTRA_BINDS")"
    TRINITY_CMD+=(--left "$R1_IN" --right "$R2_IN")
else
    TRINITY_CMD+=(--single "$R1_IN")
fi

# In-silico normalisation is Trinity's default and is usually what you want on
# deep venom-duct libraries; the toggle exists for shallow ones where it would
# discard real signal.
if [ "$NO_NORMALIZE" = "true" ]; then
    TRINITY_CMD+=(--no_normalize_reads)
fi

echo "[trinity] cmd: ${TRINITY_CMD[*]}"

t_start=$(date +%s)
# Trinity's trailing helper steps have historically exited non-zero *after*
# writing a correct assembly. Capture the status instead of letting `set -e`
# abort, then decide from the artifact itself.
set +e
bash "${HERE}/container_run.sh" "$RUNTIME" "$IMAGE" "$PLATFORM" \
     "$CACHE_DIR" "$EXTRA_BINDS" -- "${TRINITY_CMD[@]}"
RC=$?
set -e
t_end=$(date +%s)
ELAPSED=$(( t_end - t_start ))
echo "[trinity] container exit=$RC elapsed=${ELAPSED}s"

# ---- Distinguish environment faults from assembly failures -----------------
# container_run.sh reserves 10 (no runtime), 11 (daemon unreachable) and 12
# (image pull failed) for environment problems. Reporting those as the generic
# "no assembly produced" (exit 3) sends the operator hunting through Trinity's
# log for a biological explanation of what is really a stopped Docker Desktop.
# Same rationale as in 06_run_transrate.sh.
case "$RC" in
    10|11|12)
        echo "[trinity] ERROR: container runtime unavailable (exit=$RC)." >&2
        echo "[trinity] macOS: start Docker Desktop and verify with" \
             "\`docker run --rm $IMAGE Trinity --version\`." >&2
        echo "[trinity] HPC: load the apptainer/singularity module." >&2
        exit "$RC"
        ;;
esac

if [ ! -s "$NATIVE_OUT_HOST" ]; then
    NATIVE_OUT_HOST=$( { find "$OUTDIR_HOST" -maxdepth 2 -type f \
                          -name '*Trinity*.fasta' -print -quit; } 2>/dev/null || true )
fi

if [ -s "${NATIVE_OUT_HOST:-}" ]; then
    [ "$RC" -ne 0 ] && echo "[trinity] WARN: exit=$RC but assembly present;" \
                           "treating as success"
else
    echo "[trinity] ERROR: no assembly produced (exit=$RC)" >&2
    exit 3
fi

mv "$NATIVE_OUT_HOST" "$CONTIGS_OUT"
rm -f "${OUTDIR_HOST}"/*.gene_trans_map 2>/dev/null || true
rm -rf "${OUTDIR_HOST}/trinity_out_dir" 2>/dev/null || true

read -r N_CONTIGS N_BASES <<< "$(fasta_stats "$CONTIGS_OUT")"
{
    printf 'metric\tvalue\n'
    printf 'assembler\ttrinity\n'
    printf 'assembler_version\t%s\n' "${IMAGE##*:}"
    printf 'container_image\t%s\n'   "$IMAGE"
    printf 'sample_id\t%s\n'         "$SAMPLE_ID"
    printf 'layout\t%s\n'            "$LAYOUT"
    printf 'threads\t%s\n'           "$THREADS"
    printf 'max_memory_gb\t%s\n'     "$MEM_GB"
    printf 'no_normalize_reads\t%s\n' "$NO_NORMALIZE"
    printf 'elapsed_seconds\t%s\n'   "$ELAPSED"
    printf 'n_contigs\t%s\n'         "$N_CONTIGS"
    printf 'n_bases\t%s\n'           "$N_BASES"
    printf 'contigs_fasta\t%s\n'     "$CONTIGS_OUT"
} > "$SUMMARY"

echo "[trinity] done: contigs=$N_CONTIGS bases=$N_BASES"
