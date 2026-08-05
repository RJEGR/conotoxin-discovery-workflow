#!/usr/bin/env bash
# =============================================================================
# 23_stage7_predict_orfs.sh  -  Stage 7 R2: TransDecoder ORF prediction with
#                                DIAMOND-supported retention.
# =============================================================================
# Runs the canonical 3-step TransDecoder pipeline for one assembly:
#   1. TransDecoder.LongOrfs   -m <min_orf_length>
#   2. DIAMOND blastp of longest_orfs.pep vs UniProt Conidae (support hits)
#   3. TransDecoder.Predict    --retain_blastp_hits [--single_best_only]
#
# INPUTS (positional):
#   $1  contigs_fa        : /path/to/good.transrate.fa (from Stage 4)
#   $2  db_dmnd           : /path/to/uniprot.dmnd     (from resources snapshot)
#   $3  outdir            : /path/to/predict/          (created if missing)
#   $4  min_orf_length    : integer aa (config.stage7.predict.min_orf_length)
#   $5  support_evalue    : e-value cap for support blast (e.g. 1e-5)
#   $6  max_target_seqs   : diamond --max-target-seqs (1 recommended)
#   $7  threads           : integer
#   $8  single_best_only  : "true" | "false"
#   $9  no_refine_starts  : "true" | "false"  (see NOTES)
#
# OUTPUTS (written into $outdir):
#   transdecoder.pep     : predicted peptides (canonical name, downstream input)
#   transdecoder.cds     : CDS sequences
#   transdecoder.gff3    : coords on original contigs
#   transdecoder.bed     : BED coords
#   blastp.support.outfmt6 : full DIAMOND output (audit)
#   predict.summary.tsv    : n_orfs, n_blast_hits, params (audit)
#
# NOTES
# - TransDecoder writes its outputs with the input basename prefix in the CWD,
#   NOT in --output_dir. We `cd $outdir` and rename after the run.
# - `--more-sensitive` on DIAMOND (or `--ultra-sensitive` for slower/better);
#   `--sensitive` is the default modest bump. Using `--more-sensitive` here
#   because conotoxin ORFs are short and low-complexity, where more
#   sensitivity pays off substantially.
# =============================================================================

set -euo pipefail

CONTIGS_FA="$1"
DB_DMND="$2"
OUTDIR="$3"
MIN_ORF="$4"
SUPPORT_EVALUE="$5"
MAX_TARGET_SEQS="$6"
THREADS="$7"
SINGLE_BEST="$8"
NO_REFINE_STARTS="${9:-true}"

echo "[predict_orfs] contigs=$CONTIGS_FA outdir=$OUTDIR min_orf=$MIN_ORF"

mkdir -p "$OUTDIR"

# Absolute path — we cd into $OUTDIR to keep TransDecoder's cwd-relative
# side outputs contained (see NOTES above).
CONTIGS_ABS=$(cd "$(dirname "$CONTIGS_FA")" && pwd)/$(basename "$CONTIGS_FA")
DB_ABS=$(cd "$(dirname "$DB_DMND")" && pwd)/$(basename "$DB_DMND")

cd "$OUTDIR"

# Copy contigs locally with a canonical name; TransDecoder writes outputs
# prefixed with this basename.
cp "$CONTIGS_ABS" input.contigs.fasta

# ---- Step 1: TransDecoder.LongOrfs ---------------------------------------
TransDecoder.LongOrfs -t input.contigs.fasta -m "$MIN_ORF" --output_dir . \
    > longorfs.log 2>&1

# LongOrfs emits <basename>.transdecoder_dir/longest_orfs.pep
LONGEST_PEP="input.contigs.fasta.transdecoder_dir/longest_orfs.pep"
if [ ! -s "$LONGEST_PEP" ]; then
    echo "[predict_orfs] WARN: no long ORFs from LongOrfs; emitting empty outputs" >&2
    : > transdecoder.pep
    : > transdecoder.cds
    : > transdecoder.gff3
    : > transdecoder.bed
    : > blastp.support.outfmt6
    printf 'metric\tvalue\nn_orfs\t0\nn_blast_hits\t0\nmin_orf_length\t%s\nsingle_best_only\t%s\n' \
        "$MIN_ORF" "$SINGLE_BEST" > predict.summary.tsv
    rm -rf input.contigs.fasta.transdecoder_dir* input.contigs.fasta pipeliner.*.cmds 2>/dev/null || true
    exit 0
fi

# ---- Step 2: DIAMOND blastp support --------------------------------------
diamond blastp \
    --query "$LONGEST_PEP" \
    --db "$DB_ABS" \
    --outfmt 6 \
    --evalue "$SUPPORT_EVALUE" \
    --max-target-seqs "$MAX_TARGET_SEQS" \
    --threads "$THREADS" \
    --more-sensitive \
    --out blastp.support.outfmt6 \
    > diamond.log 2>&1

# ---- Step 3: TransDecoder.Predict ----------------------------------------
# `--no_refine_starts` skips the PWM-based start-codon refinement. That
# refinement trains a PWM on the input's own long ORFs and dies with
#   "Error, len(target_sequence)=NN and pwm length = 0"  (PerlLib/PWM.pm:337)
# whenever there are too few ORFs to train on. It is applied uniformly across
# every assembly (config stage7.predict.no_refine_starts) so ORF sets stay
# comparable across folds and assemblers; see the config comment for the
# full rationale.
predict_args=(-t input.contigs.fasta --retain_blastp_hits blastp.support.outfmt6 --output_dir .)
if [ "$SINGLE_BEST" = "true" ]; then
    predict_args+=(--single_best_only)
fi
if [ "$NO_REFINE_STARTS" = "true" ]; then
    predict_args+=(--no_refine_starts)
fi
TransDecoder.Predict "${predict_args[@]}" > predict.log 2>&1

# ---- Canonicalize output names -------------------------------------------
BASENAME="input.contigs.fasta"
for ext in pep cds gff3 bed; do
    if [ -e "${BASENAME}.transdecoder.${ext}" ]; then
        mv "${BASENAME}.transdecoder.${ext}" "transdecoder.${ext}"
    else
        # Predict emits empty outputs if nothing survives filtering
        : > "transdecoder.${ext}"
    fi
done

# ---- Summary -------------------------------------------------------------
N_ORFS=$(grep -c "^>" transdecoder.pep || echo 0)
N_HITS=$(wc -l < blastp.support.outfmt6 | tr -d ' ')
{
    printf 'metric\tvalue\n'
    printf 'n_orfs\t%s\n' "$N_ORFS"
    printf 'n_blast_hits\t%s\n' "$N_HITS"
    printf 'min_orf_length\t%s\n' "$MIN_ORF"
    printf 'single_best_only\t%s\n' "$SINGLE_BEST"
    printf 'no_refine_starts\t%s\n' "$NO_REFINE_STARTS"
    printf 'support_evalue\t%s\n' "$SUPPORT_EVALUE"
    printf 'max_target_seqs\t%s\n' "$MAX_TARGET_SEQS"
} > predict.summary.tsv

# ---- Cleanup -------------------------------------------------------------
# Keep the *.log files for debugging; drop the local contigs copy and
# TransDecoder's bulky intermediate checkpoint dir.
rm -rf input.contigs.fasta \
       input.contigs.fasta.transdecoder_dir \
       input.contigs.fasta.transdecoder_dir.__checkpoints \
       input.contigs.fasta.transdecoder_dir.__checkpoints_longorfs \
       pipeliner.*.cmds 2>/dev/null || true

echo "[predict_orfs] done. n_orfs=$N_ORFS n_blast_hits=$N_HITS"
