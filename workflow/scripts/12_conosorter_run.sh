#!/usr/bin/env bash
# =============================================================================
# 25_stage7_conosorter_run.sh  -  Stage 7 R4a: run ConoSorter v1.1
#                                  (RegEx + pHMM in a single invocation)
# =============================================================================
# ConoSorter's `ConoSorter -p <input.fa>` runs BOTH detection modes in one go
# and emits `<basename>_Regex.tab` and `<basename>_pHMM.tab` in the current
# working directory. It also writes many `.temp` intermediates, so we run each
# fold in its own tempdir to avoid concurrency collisions.
#
# INPUTS (positional):
#   $1  nr_pep_fasta   : /path/to/dedup/nr.pep.fasta       (from R3)
#   $2  cs_zip         : /path/to/ConoSorter_v1.1.zip      (from user drop)
#   $3  outdir         : /path/to/conosorter/              (created if missing)
#
# OUTPUTS (written into $outdir):
#   Regex.tab          : ConoSorter RegEx-mode table
#   pHMM.tab           : ConoSorter pHMM-mode table
#   conosorter.log     : combined stdout+stderr
#   run.summary.tsv    : n_regex_hits, n_phmm_hits, wall_seconds
#
# NOTES
# - Lazy unzip: if `resources/software/conosorter/ConoSorter_v1.1/ConoSorter`
#   does not exist yet, we `unzip -q -o` next to the ZIP. Idempotent.
# - Per-fold tempdir: `mktemp -d` under the outdir; we `cp -R` the entire
#   ConoSorter install into it. Cost: ~7 MB copy per fold. Alternative
#   (symlinks) breaks ConoSorter's `./transl` resolution on some macOS
#   setups; the copy is cheap and robust.
# - Input FASTA linearization: `seqkit seq -w 0` because ConoSorter's `transl`
#   awk parser expects single-line sequences (line-based FS).
# =============================================================================

set -euo pipefail

NR_PEP="$1"
CS_ZIP="$2"
OUTDIR="$3"

echo "[conosorter_run] nr_pep=$NR_PEP outdir=$OUTDIR"

# ---- Resolve absolute paths BEFORE cd ------------------------------------
NR_PEP_ABS=$(cd "$(dirname "$NR_PEP")" && pwd)/$(basename "$NR_PEP")
CS_ZIP_ABS=$(cd "$(dirname "$CS_ZIP")" && pwd)/$(basename "$CS_ZIP")
CS_INSTALL_DIR="$(dirname "$CS_ZIP_ABS")/ConoSorter_v1.1"

# ---- Lazy unzip (idempotent) ---------------------------------------------
if [ ! -x "$CS_INSTALL_DIR/ConoSorter" ]; then
    echo "[conosorter_run] Unzipping ConoSorter v1.1 (first run)..."
    (cd "$(dirname "$CS_ZIP_ABS")" && unzip -q -o "$CS_ZIP_ABS")
    # Ensure executable bits (some hosts strip them on unzip)
    chmod +x "$CS_INSTALL_DIR/ConoSorter" \
             "$CS_INSTALL_DIR/transl" \
             "$CS_INSTALL_DIR/regex_p" \
             "$CS_INSTALL_DIR/regex_d" \
             "$CS_INSTALL_DIR/pHMM" \
             "$CS_INSTALL_DIR/layout"
    # Drop the macOS metadata folder if unzip pulled it out
    rm -rf "$(dirname "$CS_ZIP_ABS")/__MACOSX" 2>/dev/null || true
fi

# Sanity: essentials must be in place
for f in ConoSorter transl regex_p pHMM layout models/HMM_superfam_db.h3m; do
    if [ ! -e "$CS_INSTALL_DIR/$f" ]; then
        echo "[conosorter_run] ERROR: missing $CS_INSTALL_DIR/$f after unzip" >&2
        exit 2
    fi
done

# ---- Per-fold tempdir with a private copy of ConoSorter ------------------
mkdir -p "$OUTDIR"
TMPDIR_LOCAL=$(mktemp -d "$OUTDIR/cs_XXXXXX")
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT

# Copy install (~7 MB) into the tempdir
cp -R "$CS_INSTALL_DIR/." "$TMPDIR_LOCAL/"

# ---- Patch hardcoded shebangs (2026-07-24 fix) ---------------------------
# ConoSorter's gawk scripts (transl, regex_p, regex_d, pHMM, layout) ship
# with `#! /usr/local/bin/gawk -f` which is not the gawk shipped by the
# conda env (it lives under .snakemake/conda/.../bin/gawk). Rewrite the
# shebang to `#!/usr/bin/env -S gawk -f` inside the per-fold tempdir so the
# install directory stays untouched and every invocation re-patches from a
# clean copy. `env -S` (split-string) has been in macOS + GNU coreutils
# since 2019 and is required to pass the `-f` argument through `env`.
for f in transl regex_p regex_d pHMM layout; do
    if [ -e "$TMPDIR_LOCAL/$f" ]; then
        # In-place edit; safe because $TMPDIR_LOCAL is unique per invocation
        # macOS `sed -i ''` vs GNU `sed -i` differ, so use portable form
        awk 'NR==1 { sub(/^#! ?\/usr\/local\/bin\/gawk -f/, "#!/usr/bin/env -S gawk -f") } { print }' \
            "$TMPDIR_LOCAL/$f" > "$TMPDIR_LOCAL/$f.patched"
        mv "$TMPDIR_LOCAL/$f.patched" "$TMPDIR_LOCAL/$f"
        chmod +x "$TMPDIR_LOCAL/$f"
    fi
done

# Linearize input FASTA into the tempdir (ConoSorter's transl expects single-line seqs)
seqkit seq -w 0 "$NR_PEP_ABS" > "$TMPDIR_LOCAL/input.fasta"

# ---- Run ConoSorter ------------------------------------------------------
LOG="$OUTDIR/conosorter.log"
START=$(date +%s)
if ! (cd "$TMPDIR_LOCAL" && bash ./ConoSorter -p input.fasta) > "$LOG" 2>&1; then
    echo "[conosorter_run] ConoSorter exited non-zero — see $LOG" >&2
    exit 3
fi
END=$(date +%s)
ELAPSED=$((END - START))

# ---- Move canonicalised outputs into $OUTDIR -----------------------------
# ConoSorter renames Regex.tab -> <basename>_Regex.tab; basename here is "input"
if [ -e "$TMPDIR_LOCAL/input_Regex.tab" ]; then
    mv "$TMPDIR_LOCAL/input_Regex.tab" "$OUTDIR/Regex.tab"
else
    : > "$OUTDIR/Regex.tab"
fi
if [ -e "$TMPDIR_LOCAL/input_pHMM.tab" ]; then
    mv "$TMPDIR_LOCAL/input_pHMM.tab" "$OUTDIR/pHMM.tab"
else
    : > "$OUTDIR/pHMM.tab"
fi

# ---- Summary -------------------------------------------------------------
# Header row is present in both tabs, so subtract 1 for the data-row count
N_REGEX=$(wc -l < "$OUTDIR/Regex.tab" | tr -d ' ')
N_REGEX=$(( N_REGEX > 0 ? N_REGEX - 1 : 0 ))
N_PHMM=$(wc -l < "$OUTDIR/pHMM.tab" | tr -d ' ')
N_PHMM=$(( N_PHMM > 0 ? N_PHMM - 1 : 0 ))
{
    printf 'metric\tvalue\n'
    printf 'n_regex_hits\t%s\n' "$N_REGEX"
    printf 'n_phmm_hits\t%s\n'  "$N_PHMM"
    printf 'wall_seconds\t%s\n' "$ELAPSED"
    printf 'conosorter_version\t1.1\n'
} > "$OUTDIR/run.summary.tsv"

echo "[conosorter_run] done. n_regex=$N_REGEX n_phmm=$N_PHMM wall=${ELAPSED}s"
