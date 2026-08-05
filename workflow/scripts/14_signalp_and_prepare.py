#!/usr/bin/env python3
# =============================================================================
# 14_signalp_and_prepare.py  -  Stage 02 R5: SignalP 6 + candidate preparation
# =============================================================================
# Adapted from the benchmark's 27_stage7_signalp_and_prepare.py. The filtering,
# containment and split logic is unchanged; what differs is HOW SignalP is
# obtained, because this pipeline must also run where the DTU academic tarball
# is not available.
#
# TWO EXECUTION ROUTES, tried in order
#
#   1. LOCAL (preferred). A DTU tarball matching `--signalp-tarball` is
#      extracted and pip-installed into the active conda env, then `signalp6`
#      runs in-process. Advantages: the exact DTU release is pinned, there is
#      no container overhead per invocation, and the model weights stay warm
#      between samples.
#
#   2. CONTAINER (fallback). `--signalp-container` is executed through
#      container_run.sh. Used automatically when no tarball is present, so a
#      collaborator can run the pipeline without first clearing the DTU
#      licence. The trade-off is real and recorded in the audit: the image is
#      third-party, its SignalP build is not version-pinned by us, and results
#      may differ from a local install if the image is ever rebuilt.
#
# Whichever route runs is written to the audit as `signalp_route`, so an
# archived result always states how its signal peptides were called.
#
# TARBALL LAYOUT (verified on 6.0i). The DTU archive is NESTED: the
# pip-installable directory is at `signalp6_fast/signalp-6-package/`, not the
# archive root, and the model weights ship separately under that directory's
# `models/` and must be copied to `$SIGNALP_DIR/model_weights/` after install.
# numpy must be held below 2.0 because SignalP pins torch<2, which was built
# against the numpy 1.x ABI and aborts with "Numpy is not available" otherwise.
#
# OUTPUT FORMAT. SignalP 6 renamed 5.0b's `--format short` to `--format txt`,
# and its `Prediction` column now carries short tags (`SP`, `OTHER`, `LIPO`,
# `TAT`) rather than the full `SP(Sec/SPI)` label — `SP(Sec/SPI)` is now the
# name of a probability column. The parser accepts both spellings.
# =============================================================================
# NOTE: no `from __future__ import annotations` — see project gotchas.
# =============================================================================

import argparse
import csv
import glob
import re
import shutil
import subprocess
import sys
from pathlib import Path

import openpyxl

SEP = "|"
SIGNALP_POSITIVE = {"SP", "SP(Sec/SPI)"}

SHARED_COLS = ["Protein Sequence", "# A.A", "# Cysteine(s)",
               "% Hydrophobicity (Signal)", "Sequence Frequency"]

PRINCIPAL_ORDER = [
    "seq_id", "Protein Sequence", "# A.A", "# Cysteine(s)",
    "% Hydrophobicity (Signal)", "Sequence Frequency",
    "regex__Superfamily (Signal)", "phmm__Superfamily (Signal)",
    "regex__Conopeptide", "phmm__Conopeptide",
    "regex__Read Name", "phmm__Read Name",
]

EMPTIES = {"", "-", "na", "n/a", "none", "null", "nan"}


# ---------------------------------------------------------------------------
# Route 1: local pip install
# ---------------------------------------------------------------------------
def try_install_local(tarball_glob: str):
    """Return the signalp6 executable path, or None if no tarball is present."""
    exe = shutil.which("signalp6")
    if exe:
        return exe
    matches = sorted(glob.glob(tarball_glob)) if tarball_glob else []
    if not matches:
        return None
    if len(matches) > 1:
        raise SystemExit(f"multiple SignalP tarballs match {tarball_glob}: "
                         f"{matches}; pin an exact filename")
    tar = matches[0]
    tar_dir = Path(tar).parent
    extract_root = tar_dir / "_extracted"
    extract_root.mkdir(exist_ok=True)
    print(f"[signalp] extracting {tar}", file=sys.stderr)
    subprocess.run(["tar", "-xzf", tar, "-C", str(extract_root)], check=True)

    pkgs = list(extract_root.glob("*/signalp-6-package"))
    if not pkgs:
        raise SystemExit(f"unexpected tarball layout under {extract_root}")
    pkg = pkgs[0]

    print(f"[signalp] pip install {pkg}", file=sys.stderr)
    subprocess.run([sys.executable, "-m", "pip", "install",
                    "--no-cache-dir", str(pkg)], check=True)
    # torch<2 was built against the numpy 1.x ABI; a transitive upgrade to
    # numpy>=2 makes every `.numpy()` call abort at inference time.
    subprocess.run([sys.executable, "-m", "pip", "install", "--no-cache-dir",
                    "--force-reinstall", "--no-deps", "numpy<2"], check=True)

    r = subprocess.run(
        [sys.executable, "-c",
         "import signalp, os; print(os.path.dirname(signalp.__file__))"],
        check=True, capture_output=True, text=True)
    weights_dst = Path(r.stdout.strip()) / "model_weights"
    weights_dst.mkdir(exist_ok=True)
    models_src = pkg / "models"
    if models_src.is_dir():
        n = 0
        for item in models_src.iterdir():
            dst = weights_dst / item.name
            if item.is_file() and not dst.exists():
                shutil.copy2(item, dst)
                n += 1
        print(f"[signalp] copied {n} weight file(s)", file=sys.stderr)

    exe = shutil.which("signalp6")
    if not exe:
        raise SystemExit("pip install succeeded but signalp6 is not on PATH")
    return exe


# ---------------------------------------------------------------------------
# SignalP invocation (either route)
# ---------------------------------------------------------------------------
def run_signalp(fasta: Path, outdir: Path, organism: str, mode: str,
                fmt: str, exe, container_args):
    outdir.mkdir(parents=True, exist_ok=True)
    existing = outdir / "prediction_results.txt"
    if existing.exists() and existing.stat().st_size:
        print(f"[signalp] reusing {existing} (delete to force re-inference)",
              file=sys.stderr)
        return existing, "reused"

    base = ["--fastafile", str(fasta), "--organism", organism,
            "--output_dir", str(outdir), "--format", fmt, "--mode", mode]

    if exe:
        cmd = [exe] + base
        route = "local"
    else:
        runtime, image, platform, cache, binds, helper = container_args
        # container_run.sh binds the workdir to /work, so workdir-relative
        # paths become /work/<path> inside the container.
        cmd = ["bash", helper, runtime, image, platform, cache, binds, "--",
               "signalp6",
               "--fastafile", f"/work/{fasta}",
               "--organism", organism,
               "--output_dir", f"/work/{outdir}",
               "--format", fmt, "--mode", mode]
        route = "container"

    print(f"[signalp] route={route}: {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd, check=True)

    pred = outdir / "prediction_results.txt"
    if not pred.exists():
        cands = list(outdir.glob("*results*.txt"))
        if not cands:
            raise SystemExit(f"SignalP finished but no results file in {outdir}")
        pred = cands[0]
    return pred, route


def parse_positive_ids(path: Path) -> set:
    pos, header_found, id_i, pred_i = set(), False, None, None
    with path.open(encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n\r")
            if not line:
                continue
            if not header_found and line.startswith("# ID") and "Prediction" in line:
                cols = line.split("\t")
                cols[0] = re.sub(r"^\#\s*", "", cols[0])
                id_i, pred_i = cols.index("ID"), cols.index("Prediction")
                header_found = True
                continue
            if line.startswith("#") or not header_found:
                continue
            parts = line.split("\t")
            if len(parts) > max(id_i, pred_i) and parts[pred_i].strip() in SIGNALP_POSITIVE:
                pos.add(parts[id_i].strip())
    return pos


# ---------------------------------------------------------------------------
# Table handling (unchanged from the benchmark)
# ---------------------------------------------------------------------------
def read_merge_tab(path: Path):
    with path.open(encoding="utf-8", errors="replace") as fh:
        headers = fh.readline().rstrip("\n\r").split(SEP)
        rows = []
        for raw in fh:
            line = raw.rstrip("\n\r")
            if not line:
                continue
            parts = line.split(SEP)
            if len(parts) < len(headers):
                parts += [""] * (len(headers) - len(parts))
            elif len(parts) > len(headers):
                parts = parts[:len(headers) - 1] + [SEP.join(parts[len(headers) - 1:])]
            rows.append({h: p.strip() for h, p in zip(headers, parts)})
    return headers, rows


def write_pipe(path: Path, headers, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write(SEP.join(headers) + "\n")
        for r in rows:
            fh.write(SEP.join(str(r.get(h, "")).replace(SEP, "_")
                              for h in headers) + "\n")


def containment_filter(rows, seq_col="Protein Sequence"):
    """Drop sequences that are exact substrings of a longer retained one."""
    kept, kept_seqs = [], []
    for r in sorted(rows, key=lambda x: len(x.get(seq_col, "")), reverse=True):
        s = r.get(seq_col, "").upper()
        if not s or any(s in longer for longer in kept_seqs):
            continue
        kept.append(r)
        kept_seqs.append(s)
    return sorted(kept, key=lambda r: r.get("seq_id", ""))


def is_assigned(row) -> bool:
    for k in ("regex__Superfamily (Signal)", "phmm__Superfamily (Signal)"):
        v = str(row.get(k, "")).strip().lower()
        if v and v not in EMPTIES:
            return True
    return False


def write_xlsx(headers, rows, out_xlsx: Path):
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "principal"
    principal = [c for c in PRINCIPAL_ORDER if c in headers]
    movidas = [c for c in headers if c not in principal]
    ws.append(principal)
    for r in rows:
        ws.append([r.get(c, "") for c in principal])
    ws2 = wb.create_sheet("movidas")
    ws2.append(["seq_id"] + movidas)
    for r in rows:
        ws2.append([r.get("seq_id", "")] + [r.get(c, "") for c in movidas])
    wb.save(out_xlsx)


def write_fasta(rows, out_fa: Path):
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with out_fa.open("w") as fh:
        for r in rows:
            sid = r.get("seq_id", "").strip()
            seq = r.get("Protein Sequence", "").strip().upper().replace("*", "")
            if sid and seq:
                fh.write(f">{sid}\n{seq}\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--merge-tab",   required=True, type=Path)
    ap.add_argument("--merge-fasta", required=True, type=Path)
    ap.add_argument("--outdir",      required=True, type=Path)
    ap.add_argument("--signalp-tarball",   default="")
    ap.add_argument("--signalp-container", default="")
    ap.add_argument("--organism", default="eukarya")
    ap.add_argument("--mode",     default="fast")
    ap.add_argument("--format",   default="txt")
    ap.add_argument("--container-runtime",  default="auto")
    ap.add_argument("--container-platform", default="")
    ap.add_argument("--container-cache",    default="resources/containers")
    ap.add_argument("--container-binds",    default="")
    ap.add_argument("--container-helper",
                    default="workflow/scripts/container_run.sh")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    exe = try_install_local(args.signalp_tarball)
    if exe is None and not args.signalp_container:
        raise SystemExit(
            "no SignalP source: neither a DTU tarball at "
            f"'{args.signalp_tarball}' nor a container image configured")

    pred, route = run_signalp(
        args.merge_fasta, args.outdir / "signalp_output",
        args.organism, args.mode, args.format, exe,
        (args.container_runtime, args.signalp_container,
         args.container_platform, args.container_cache,
         args.container_binds, args.container_helper))

    pos_ids = parse_positive_ids(pred)
    headers, merge_rows = read_merge_tab(args.merge_tab)
    filtered = [r for r in merge_rows if r.get("seq_id") in pos_ids]
    write_pipe(args.outdir / "filtered.tab", headers, filtered)

    kept = containment_filter(filtered)
    assigned = [r for r in kept if is_assigned(r)]
    unassigned = [r for r in kept if not is_assigned(r)]

    write_xlsx(headers, kept,       args.outdir / "reorganized.xlsx")
    write_xlsx(headers, assigned,   args.outdir / "positivo_superfamilia.xlsx")
    write_xlsx(headers, unassigned, args.outdir / "negativo_superfamilia.xlsx")
    write_fasta(assigned,   args.outdir / "positivo.fasta")
    write_fasta(unassigned, args.outdir / "negativo.fasta")

    audit = args.outdir / "signalp_prepare.audit.tsv"
    with audit.open("w", encoding="utf-8") as fh:
        fh.write("metric\tvalue\n")
        for k, v in (("signalp_route", route),
                     ("signalp_organism", args.organism),
                     ("signalp_mode", args.mode),
                     ("n_merge_rows", len(merge_rows)),
                     ("n_signalp_positive", len(pos_ids)),
                     ("n_filtered", len(filtered)),
                     ("n_after_containment", len(kept)),
                     ("n_assigned", len(assigned)),
                     ("n_unassigned", len(unassigned))):
            fh.write(f"{k}\t{v}\n")

    print(f"[signalp] route={route} merge={len(merge_rows)} "
          f"sp_pos={len(pos_ids)} kept={len(kept)} "
          f"assigned={len(assigned)} unassigned={len(unassigned)}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
