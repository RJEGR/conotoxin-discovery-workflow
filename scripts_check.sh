#!/usr/bin/env bash
# =============================================================================
# scripts_check.sh  -  static validation of the pipeline, no data required
# =============================================================================
# Run this after any edit, and on a fresh clone before the first real run. It
# catches the class of errors that would otherwise surface an hour into a job:
# a syntax slip, a missing vendored script, a malformed YAML, a sample sheet
# indented with spaces.
#
#   bash scripts_check.sh
#
# Exits non-zero if anything fails, so it can gate a CI job.
# =============================================================================

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

FAILED=0
pass() { printf '  \033[32mOK\033[0m   %s\n' "$1"; }
fail() { printf '  \033[31mFAIL\033[0m %s\n' "$1"; FAILED=1; }

echo "== bash syntax =="
for f in workflow/scripts/*.sh; do
    [ -e "$f" ] || continue
    if bash -n "$f" 2>/dev/null; then pass "$f"; else fail "$f"; fi
done

echo
echo "== python syntax =="
for f in workflow/scripts/*.py workflow/scripts/lib/*.py; do
    [ -e "$f" ] || continue
    if python3 -m py_compile "$f" 2>/dev/null; then pass "$f"; else fail "$f"; fi
done

echo
echo "== yaml syntax =="
for f in config/config.yaml workflow/envs/*.yaml profiles/*/config.yaml; do
    [ -e "$f" ] || continue
    if python3 -c "import yaml,sys; yaml.safe_load(open(sys.argv[1]))" "$f" 2>/dev/null; then
        pass "$f"
    else
        fail "$f"
    fi
done

echo
echo "== vendored scripts present =="
# These are copied from ../workflow/scripts by the setup step; a missing one
# means the copy was skipped and the DAG will fail at that rule.
for f in \
    workflow/scripts/lib/legacy11_normalize.py \
    workflow/scripts/lib/legacy12_classify.py \
    workflow/scripts/10_predict_orfs.sh \
    workflow/scripts/11_dedup_merge_pep.sh \
    workflow/scripts/12_conosorter_run.sh \
    workflow/scripts/13_merge_conosorter.py \
    workflow/scripts/15_diamond_annotate.sh \
    workflow/scripts/16_parse_besthits.py \
    workflow/scripts/17_annotate_labels.py \
    workflow/scripts/18_normalize_final.py \
    workflow/scripts/19_classify_confidence.py ; do
    if [ -s "$f" ]; then pass "$f"; else fail "$f  (missing — re-run the cp step in the README)"; fi
done

echo
echo "== sample sheet =="
if [ -s config/samples.tsv ]; then
    # A sheet indented with spaces parses as a single column and produces a
    # confusing "missing required column" failure later.
    if head -1 config/samples.tsv | grep -q $'\t'; then
        pass "config/samples.tsv is tab-separated"
    else
        fail "config/samples.tsv header has no tabs (check with: cat -A config/samples.tsv)"
    fi
    ncol=$(head -1 config/samples.tsv | awk -F'\t' '{print NF}')
    nbad=$(awk -F'\t' -v n="$ncol" 'NR>1 && NF>0 && NF>n {c++} END{print c+0}' config/samples.tsv)
    if [ "$nbad" -eq 0 ]; then
        pass "all rows have <= $ncol fields"
    else
        fail "$nbad row(s) have more fields than the header"
    fi
else
    printf '  \033[33mSKIP\033[0m config/samples.tsv not created yet\n'
fi

echo
echo "== container runtime =="
if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
    mem=$(docker info 2>/dev/null | awk -F': *' '/Total Memory/{print $2}')
    pass "docker reachable (memory: ${mem:-unknown})"
    printf '       keep assembly.assemblers.*.max_memory_gb at ~75%% of that\n'
elif command -v apptainer >/dev/null 2>&1; then
    pass "apptainer available"
elif command -v singularity >/dev/null 2>&1; then
    pass "singularity available"
else
    printf '  \033[33mWARN\033[0m no container runtime: Trinity, TransRate and the\n'
    printf '       SignalP fallback will fail. Start Docker Desktop or load the\n'
    printf '       apptainer module.\n'
fi

echo
if [ "$FAILED" -eq 0 ]; then
    printf '\033[32mAll static checks passed.\033[0m\n'
else
    printf '\033[31mSome checks failed (see above).\033[0m\n'
fi
exit "$FAILED"
