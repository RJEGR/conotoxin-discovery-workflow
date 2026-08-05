#!/usr/bin/env bash
# =============================================================================
# container_run.sh  -  portable container invocation (Docker | Apptainer)
# =============================================================================
# Single entry point used by every containerised rule (Trinity, TransRate,
# SignalP). It hides the two runtimes behind one calling convention so the
# same rule code runs on a macOS laptop with Docker Desktop and on an HPC node
# where only Apptainer/Singularity exists (Docker needs root and is almost
# never available on shared clusters).
#
# USAGE
#   container_run.sh <runtime> <image> <platform> <cache_dir> <extra_binds> -- CMD...
#
#     runtime      auto | docker | apptainer | singularity
#     image        docker-style reference, e.g. trinityrnaseq/trinityrnaseq:2.15.2
#     platform     value for `docker run --platform` (ignored by apptainer and
#                  on native linux/amd64 hosts); pass "" to omit
#     cache_dir    where .sif images are cached (apptainer only)
#     extra_binds  comma-separated host:container pairs, or "" for none
#     CMD...       the command to execute INSIDE the container
#
# PATH MODEL (identical under both runtimes, so rule code never branches)
#   host $PWD  <->  /work    and the container's working directory is /work.
#   Therefore any path a rule passes must be expressed RELATIVE TO $PWD and
#   prefixed with /work — use the `to_container_path` helper below.
#   Additional binds from `extra_binds` are mounted at the SAME absolute path
#   inside the container, which is what you want for read files living outside
#   the workdir (a very common HPC layout: reads on /scratch, run from /home).
#
# EXIT CODES
#   10  no usable container runtime found
#   11  docker selected but daemon unreachable
#   12  image pull failed
#   <n> otherwise: the containerised command's own exit code, passed through
#       verbatim so callers can implement tool-specific recovery.
# =============================================================================

set -euo pipefail

RUNTIME_REQ="${1:?runtime required}"
IMAGE="${2:?image required}"
PLATFORM="${3:-}"
CACHE_DIR="${4:-resources/containers}"
EXTRA_BINDS="${5:-}"
shift 5
[ "${1:-}" = "--" ] && shift

if [ "$#" -eq 0 ]; then
    echo "[container_run] ERROR: no command given after --" >&2
    exit 2
fi

# ---------------------------------------------------------------------------
# Runtime detection
# ---------------------------------------------------------------------------
detect_runtime() {
    case "$RUNTIME_REQ" in
        docker|apptainer|singularity)
            echo "$RUNTIME_REQ"; return 0 ;;
        auto) : ;;
        *)  echo "[container_run] ERROR: unknown runtime '$RUNTIME_REQ'" >&2
            exit 10 ;;
    esac
    # Preference order: docker first (fast, no .sif conversion), then the
    # rootless HPC options.
    if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
        echo docker; return 0
    fi
    if command -v apptainer >/dev/null 2>&1; then echo apptainer; return 0; fi
    if command -v singularity >/dev/null 2>&1; then echo singularity; return 0; fi
    echo "[container_run] ERROR: no usable container runtime. Install Docker" \
         "Desktop (macOS) or load an apptainer/singularity module (HPC)." >&2
    exit 10
}

RUNTIME="$(detect_runtime)"

# ---------------------------------------------------------------------------
# Native platform check: --platform emulation is only meaningful when the host
# architecture differs from the image's. Passing --platform linux/amd64 on a
# native amd64 host is harmless but noisy; on arm64 it is essential.
# ---------------------------------------------------------------------------
HOST_ARCH="$(uname -m)"
USE_PLATFORM=""
if [ -n "$PLATFORM" ] && [ "$RUNTIME" = "docker" ]; then
    case "$HOST_ARCH" in
        arm64|aarch64) USE_PLATFORM="$PLATFORM" ;;
        *)             USE_PLATFORM="" ;;
    esac
fi

HOST_WORKDIR="$(pwd)"

# ---------------------------------------------------------------------------
# Assemble bind arguments
# ---------------------------------------------------------------------------
BIND_ARGS=()
add_bind() {
    local spec="$1"
    [ -z "$spec" ] && return 0
    if [ "$RUNTIME" = "docker" ]; then
        BIND_ARGS+=(-v "$spec")
    else
        BIND_ARGS+=(--bind "$spec")
    fi
}

add_bind "${HOST_WORKDIR}:/work"
if [ -n "$EXTRA_BINDS" ]; then
    IFS=',' read -ra _binds <<< "$EXTRA_BINDS"
    for b in "${_binds[@]}"; do
        b="$(echo "$b" | xargs)"     # trim
        [ -n "$b" ] && add_bind "$b"
    done
fi

echo "[container_run] runtime=$RUNTIME image=$IMAGE host_arch=$HOST_ARCH" >&2
[ -n "$USE_PLATFORM" ] && echo "[container_run] platform=$USE_PLATFORM" >&2

# ---------------------------------------------------------------------------
# Execute
# ---------------------------------------------------------------------------
if [ "$RUNTIME" = "docker" ]; then
    if ! docker info >/dev/null 2>&1; then
        echo "[container_run] ERROR: docker daemon not reachable." >&2
        exit 11
    fi
    DOCKER_ARGS=(run --rm)
    [ -n "$USE_PLATFORM" ] && DOCKER_ARGS+=(--platform "$USE_PLATFORM")
    DOCKER_ARGS+=("${BIND_ARGS[@]}" -w /work)
    # Run as the invoking user so outputs are not root-owned on the host.
    DOCKER_ARGS+=(-u "$(id -u):$(id -g)")
    DOCKER_ARGS+=("$IMAGE" "$@")
    exec docker "${DOCKER_ARGS[@]}"
fi

# ---- Apptainer / Singularity ----------------------------------------------
# Both consume docker:// URIs and cache a .sif. We materialise the .sif
# explicitly (rather than letting the runtime re-resolve the registry on every
# invocation) so that a batch of array jobs does not hammer Docker Hub and
# trip its rate limit -- a classic HPC failure mode.
mkdir -p "$CACHE_DIR"
SIF_NAME="$(echo "$IMAGE" | tr '/:' '__').sif"
SIF_PATH="${CACHE_DIR}/${SIF_NAME}"

if [ ! -s "$SIF_PATH" ]; then
    echo "[container_run] building $SIF_PATH from docker://$IMAGE" >&2
    # Serialise concurrent builds with a lock dir: many array tasks starting at
    # once would otherwise each pull the same image and corrupt the cache.
    #
    # The lock is released EXPLICITLY rather than through `trap ... EXIT`,
    # because this script ends in `exec` — which replaces the process image
    # without firing EXIT traps, leaving the lock orphaned and stalling every
    # later invocation for the full timeout.
    LOCK="${SIF_PATH}.lock"
    # A lock older than an hour is a leftover from a killed job (SLURM
    # preemption, node failure): reclaim it rather than waiting it out.
    if [ -d "$LOCK" ] && [ -z "$(find "$LOCK" -maxdepth 0 -mmin -60 2>/dev/null)" ]; then
        echo "[container_run] removing stale lock $LOCK" >&2
        rmdir "$LOCK" 2>/dev/null || true
    fi
    for _ in $(seq 1 120); do
        if mkdir "$LOCK" 2>/dev/null; then
            build_rc=0
            if [ ! -s "$SIF_PATH" ]; then
                "$RUNTIME" build --force "$SIF_PATH" "docker://$IMAGE" >&2 || build_rc=$?
            fi
            rmdir "$LOCK" 2>/dev/null || true
            if [ "$build_rc" -ne 0 ]; then
                echo "[container_run] ERROR: failed to build $SIF_PATH" >&2
                exit 12
            fi
            break
        fi
        # Another task holds the lock; wait for it to finish.
        sleep 5
        [ -s "$SIF_PATH" ] && break
    done
fi

if [ ! -s "$SIF_PATH" ]; then
    echo "[container_run] ERROR: $SIF_PATH not available after build/wait" >&2
    exit 12
fi

APP_ARGS=(exec --cleanenv "${BIND_ARGS[@]}" --pwd /work "$SIF_PATH" "$@")
exec "$RUNTIME" "${APP_ARGS[@]}"
