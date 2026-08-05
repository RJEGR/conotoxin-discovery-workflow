#!/usr/bin/env bash
# =============================================================================
# 24_stage7_dedup_merge_pep.sh  -  Stage 7 R3: cross-assembler peptide
#                                    dereplication with MMseqs2 easy-cluster.
# =============================================================================
# Concatenates the 3 per-assembler `transdecoder.pep` outputs from R2 (with
# assembler-tagged headers to preserve provenance) and dereplicates them into
# a single non-redundant peptide FASTA (`nr.pep.fasta`) that R4 feeds to
# ConoSorter.
#
# Rationale (agreed 2026-07-24): MMseqs2 replaces the legacy CD-HIT-2D. Same
# thresholds (0.95 id, 0.90 coverage of shorter member) translated to
# MMseqs2 nomenclature (--min-seq-id + -c + --cov-mode 1).
#
# INPUTS (positional):
#   $1  outdir           : /path/to/dedup/  (created if missing)
#   $2  min_seq_id       : 0.0..1.0    (e.g. 0.95)
#   $3  coverage         : 0.0..1.0    (e.g. 0.90)
#   $4  cov_mode         : 0|1|2       (1 = coverage of target/shorter)
#   $5  sensitivity      : e.g. 7.5
#   $6  threads          : integer
#   $7..N pep_paths      : one per assembler (order does not matter; tag comes
#                           from the parent .../{asm}/predict/ dir)
#
# OUTPUTS (written into $outdir):
#   nr.pep.fasta          : dereplicated representative sequences
#   nr.cluster.tsv        : two-column {cluster_rep, member}
#   combined.pep.fasta    : provenance-tagged concatenation of inputs
#   dedup.audit.tsv       : {cluster_rep, member, source_assembler}
#   cluster_members.tsv   : histogram of members-per-cluster
#   dedup.summary.tsv     : n_input_orfs, n_nr_orfs, n_clusters, params
#   mmseqs.log            : full MMseqs2 stdout+stderr (audit)
# =============================================================================

set -euo pipefail

OUTDIR="$1"
MIN_ID="$2"
COV="$3"
COV_MODE="$4"
SENS="$5"
THREADS="$6"
shift 6
PEPS=("$@")

echo "[dedup_merge_pep] n_inputs=${#PEPS[@]} outdir=$OUTDIR min_id=$MIN_ID cov=$COV"

# ---- Resolve absolute paths BEFORE cd ------------------------------------
# Snakemake passes input paths relative to the workdir (repo root). After the
# `cd "$OUTDIR"` below those relative paths point to nothing, so we absolutize
# them here first. Same fix as in script 23.
PEPS_ABS=()
for pep in "${PEPS[@]}"; do
    if [ ! -e "$pep" ]; then
        echo "[dedup_merge_pep] WARN: input missing (pre-cd): $pep" >&2
        continue
    fi
    PEPS_ABS+=("$(cd "$(dirname "$pep")" && pwd)/$(basename "$pep")")
done

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# ---- Concatenate with assembler-tagged headers ---------------------------
# Header rewriting: prefix each seq id with `<asm>_` so cross-assembler
# provenance is recoverable from the cluster table. asm name is inferred
# from the input path: .../{asm}/predict/transdecoder.pep -> {asm}.
: > combined.pep.fasta
N_INPUT=0
for pep in "${PEPS_ABS[@]}"; do
    asm=$(basename "$(dirname "$(dirname "$pep")")")
    awk -v asm="$asm" '/^>/ { sub(/^>/, ">" asm "_"); print; next } { print }' \
        "$pep" >> combined.pep.fasta
    N_INPUT=$((N_INPUT + $(grep -c "^>" "$pep" || true)))
done

# Guard: MMseqs2 refuses to run on an empty input.
if [ ! -s combined.pep.fasta ]; then
    echo "[dedup_merge_pep] WARN: combined.pep.fasta is empty; emitting empty outputs" >&2
    : > nr.pep.fasta
    : > nr.cluster.tsv
    : > dedup.audit.tsv
    : > cluster_members.tsv
    printf 'metric\tvalue\nn_input_orfs\t0\nn_nr_orfs\t0\nn_clusters\t0\n' > dedup.summary.tsv
    exit 0
fi

# ---- MMseqs2 easy-cluster -----------------------------------------------
mmseqs easy-cluster combined.pep.fasta nr tmp_mmseqs \
    --min-seq-id "$MIN_ID" \
    -c "$COV" \
    --cov-mode "$COV_MODE" \
    -s "$SENS" \
    --threads "$THREADS" \
    --remove-tmp-files 1 \
    > mmseqs.log 2>&1

mv nr_rep_seq.fasta nr.pep.fasta
mv nr_cluster.tsv   nr.cluster.tsv

# ---- Audit tables --------------------------------------------------------
# members-per-cluster histogram
awk -F'\t' '{ c[$1]++ } END { for (k in c) print k"\t"c[k] }' nr.cluster.tsv \
    | sort -k2,2nr > cluster_members.tsv

# canonical audit TSV: {cluster_rep, member, source_assembler}
{
    printf 'cluster_rep\tmember\tsource_assembler\n'
    while IFS=$'\t' read -r rep member; do
        src=$(printf '%s' "$member" | cut -d_ -f1)
        printf '%s\t%s\t%s\n' "$rep" "$member" "$src"
    done < nr.cluster.tsv
} > dedup.audit.tsv

N_NR=$(grep -c "^>" nr.pep.fasta || echo 0)
N_CLUSTERS=$(cut -f1 nr.cluster.tsv | sort -u | wc -l | tr -d ' ')

{
    printf 'metric\tvalue\n'
    printf 'n_input_orfs\t%s\n' "$N_INPUT"
    printf 'n_nr_orfs\t%s\n' "$N_NR"
    printf 'n_clusters\t%s\n' "$N_CLUSTERS"
    printf 'min_seq_id\t%s\n' "$MIN_ID"
    printf 'coverage\t%s\n' "$COV"
    printf 'cov_mode\t%s\n' "$COV_MODE"
    printf 'sensitivity\t%s\n' "$SENS"
} > dedup.summary.tsv

# Drop the all-seqs FASTA (large; identical info in cluster.tsv + nr.pep.fasta)
rm -f nr_all_seqs.fasta 2>/dev/null || true

echo "[dedup_merge_pep] done. n_input=$N_INPUT n_nr=$N_NR n_clusters=$N_CLUSTERS"
