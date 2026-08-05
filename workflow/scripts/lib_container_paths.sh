#!/usr/bin/env bash
# =============================================================================
# lib_container_paths.sh  -  shared helpers for containerised assembly rules
# =============================================================================
# Sourced (not executed) by the Trinity and TransRate wrappers. Centralises the
# host<->container path arithmetic so the two wrappers cannot drift apart.
#
# PATH MODEL. container_run.sh always binds the host workdir to /work. A path
# therefore has two possible fates:
#
#   * INSIDE the workdir  -> rewrite as /work/<relative path>. This covers
#     everything the pipeline produces (results/, resources/) and is the
#     normal case.
#   * OUTSIDE the workdir -> keep the absolute path unchanged and require the
#     user to have bound it via `containers.extra_binds`. This is the common
#     HPC layout where raw FASTQs live on /scratch while the run happens on
#     /home, and it only arises when `preprocessing.fastp` is false (with
#     fastp on, Stage 01 consumes files under results/).
#
# `to_container_path` implements that rule and fails loudly on the second case
# when the bind is missing, because the alternative is a container that starts
# fine and then reports "file not found" from deep inside an assembler log.
# =============================================================================

# to_container_path <host_path> [extra_binds_csv]
#   Echoes the path as seen inside the container.
to_container_path() {
    local host_path="$1"
    local extra_binds="${2:-}"
    local workdir; workdir="$(pwd -P)"

    # Resolve to an absolute, symlink-free path without requiring the file to
    # exist yet (outputs are named before they are created).
    local abs
    abs="$(python -c '
import os, sys
p = os.path.abspath(sys.argv[1])
d = os.path.dirname(p)
if os.path.isdir(d):
    d = os.path.realpath(d)
print(os.path.join(d, os.path.basename(p)))
' "$host_path")"

    case "$abs" in
        "$workdir"/*)
            echo "/work/${abs#"$workdir"/}"
            return 0
            ;;
        "$workdir")
            echo "/work"
            return 0
            ;;
    esac

    # Outside the workdir: only valid if an extra bind covers it.
    if [ -n "$extra_binds" ]; then
        local IFS=','
        for b in $extra_binds; do
            local host_side="${b%%:*}"
            host_side="$(echo "$host_side" | xargs)"
            case "$abs" in
                "$host_side"/*|"$host_side")
                    # extra binds are mounted at the same absolute path
                    echo "$abs"
                    return 0
                    ;;
            esac
        done
    fi

    echo "[container_paths] ERROR: '$host_path' lies outside the workdir" \
         "($workdir) and is not covered by containers.extra_binds." \
         "Add its parent directory to extra_binds in config.yaml." >&2
    return 1
}

# decompress_if_gz <path> <scratch_dir>
#   TransRate 1.0.3 predates transparent gzip input. Echoes a path to a plain
#   FASTQ, decompressing into $scratch_dir when needed.
decompress_if_gz() {
    local src="$1" scratch="$2"
    case "$src" in
        *.gz)
            local base; base="$(basename "${src%.gz}")"
            local dst="${scratch}/${base}"
            if [ ! -s "$dst" ]; then
                if command -v pigz >/dev/null 2>&1; then
                    pigz -dc "$src" > "$dst"
                else
                    gzip -dc "$src" > "$dst"
                fi
            fi
            echo "$dst"
            ;;
        *)
            echo "$src"
            ;;
    esac
}

# fasta_stats <fasta>
#   Echoes "<n_contigs>\t<n_bases>". Used by every assembler summary.
fasta_stats() {
    local fa="$1"
    if [ ! -s "$fa" ]; then
        printf '0\t0\n'
        return 0
    fi
    awk 'BEGIN{n=0; b=0}
         /^>/ {n++; next}
         {gsub(/[[:space:]]/, ""); b += length($0)}
         END{printf "%d\t%d\n", n, b}' "$fa"
}
