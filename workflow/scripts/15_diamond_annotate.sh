#!/usr/bin/env bash
# =============================================================================
# 28_stage7_diamond_annotate.sh  -  Stage 7 R6a: DIAMOND blastp of the
#                                    SignalP-filtered candidates against the
#                                    three Conidae reference databases.
# =============================================================================
# Runs 6 searches per fold: {positivo, negativo} x {uniprot, ncbi, conoserver}.
# Parameters come verbatim from the Stage 7 <context> annotation step 3:
#   e-value 1e-3, query coverage >= 30%, --max-hsps 1, and database-specific
#   --max-target-seqs (2 for UniProt and NCBI, 3 for ConoServer).
#
# INPUTS (positional):
#   $1  pos_fasta        : signalp/positivo.fasta   (ConoSorter-assigned set)
#   $2  neg_fasta        : signalp/negativo.fasta   (ConoSorter-unassigned set)
#   $3  db_uniprot       : *.dmnd
#   $4  db_ncbi          : *.dmnd
#   $5  db_conoserver    : *.dmnd
#   $6  outdir           : annotate/blast/  (created if missing)
#   $7  evalue           : e.g. 1e-3
#   $8  min_qcov         : e.g. 30           -> --query-cover
#   $9  max_hsps         : e.g. 1
#   $10 mts_uniprot      : e.g. 2
#   $11 mts_ncbi         : e.g. 2
#   $12 mts_conoserver   : e.g. 3
#   $13 threads          : integer
#
# OUTPUTS (in $outdir):
#   {positivo,negativo}__{uniprot,ncbi,conoserver}.outfmt6   (6 files)
#   diamond.summary.tsv                                      (row counts)
#
# OUTFMT (16 columns, matches the legacy script-7 contract):
#   qseqid qlen sseqid slen pident length mismatch gapopen
#   qstart qend sstart send evalue bitscore qcovhsp scovhsp
#
# NOTES
# - An empty query FASTA is legitimate (e.g. a fold where every candidate got
#   a ConoSorter superfamily => negativo.fasta is empty). DIAMOND errors on
#   empty input, so we emit an empty .outfmt6 and continue.
# - `--more-sensitive` for the same reason as R2: conotoxin peptides are short
#   and low-complexity, where extra sensitivity materially improves recall.
# =============================================================================

set -euo pipefail

POS_FASTA="$1"
NEG_FASTA="$2"
DB_UNIPROT="$3"
DB_NCBI="$4"
DB_CONOSERVER="$5"
OUTDIR="$6"
EVALUE="$7"
MIN_QCOV="$8"
MAX_HSPS="$9"
MTS_UNIPROT="${10}"
MTS_NCBI="${11}"
MTS_CONOSERVER="${12}"
THREADS="${13}"

echo "[diamond_annotate] outdir=$OUTDIR evalue=$EVALUE qcov=$MIN_QCOV threads=$THREADS"

OUTFMT="6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp"

mkdir -p "$OUTDIR"

# ---- Resolve absolute paths BEFORE any cd (gotcha from R3) ----------------
abspath() { echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
POS_ABS=$(abspath "$POS_FASTA")
NEG_ABS=$(abspath "$NEG_FASTA")
DB_UNIPROT_ABS=$(abspath "$DB_UNIPROT")
DB_NCBI_ABS=$(abspath "$DB_NCBI")
DB_CONOSERVER_ABS=$(abspath "$DB_CONOSERVER")

#   run_one <part> <query_fasta> <db_dmnd> <max_target_seqs> <db_tag>
run_one() {
    local part="$1" query="$2" db="$3" mts="$4" tag="$5"
    local out="$OUTDIR/${part}__${tag}.outfmt6"

    if [ ! -s "$query" ]; then
        echo "[diamond_annotate] $part/$tag: empty query, emitting empty outfmt6" >&2
        : > "$out"
        return 0
    fi

    diamond blastp \
        --query "$query" \
        --db "$db" \
        --outfmt $OUTFMT \
        --evalue "$EVALUE" \
        --query-cover "$MIN_QCOV" \
        --max-hsps "$MAX_HSPS" \
        --max-target-seqs "$mts" \
        --threads "$THREADS" \
        --more-sensitive \
        --out "$out" \
        >> "$OUTDIR/diamond.log" 2>&1
}

: > "$OUTDIR/diamond.log"

run_one positivo "$POS_ABS" "$DB_UNIPROT_ABS"    "$MTS_UNIPROT"    uniprot
run_one positivo "$POS_ABS" "$DB_NCBI_ABS"       "$MTS_NCBI"       ncbi
run_one positivo "$POS_ABS" "$DB_CONOSERVER_ABS" "$MTS_CONOSERVER" conoserver
run_one negativo "$NEG_ABS" "$DB_UNIPROT_ABS"    "$MTS_UNIPROT"    uniprot
run_one negativo "$NEG_ABS" "$DB_NCBI_ABS"       "$MTS_NCBI"       ncbi
run_one negativo "$NEG_ABS" "$DB_CONOSERVER_ABS" "$MTS_CONOSERVER" conoserver

# ---- Summary -------------------------------------------------------------
{
    printf 'part\tdb\tn_hits\n'
    for part in positivo negativo; do
        for db in uniprot ncbi conoserver; do
            f="$OUTDIR/${part}__${db}.outfmt6"
            n=$(wc -l < "$f" | tr -d ' ')
            printf '%s\t%s\t%s\n' "$part" "$db" "$n"
        done
    done
} > "$OUTDIR/diamond.summary.tsv"

echo "[diamond_annotate] done."
cat "$OUTDIR/diamond.summary.tsv"
