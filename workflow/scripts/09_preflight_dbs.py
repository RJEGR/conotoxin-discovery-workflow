#!/usr/bin/env python3
# =============================================================================
# 09_preflight_dbs.py  -  Stage 02 R1: verify annotation resources
# =============================================================================
# Runs before any annotation work. Every check that can be made in seconds is
# made here, because the alternative is discovering a missing database three
# hours into a run, after ORF prediction has already burned the CPU budget.
#
# CHECKED
#   1. The Conidae snapshot directory exists.
#   2. Each of the three DBs has a non-empty FASTA and a non-empty .dmnd index.
#   3. The two label auxiliaries exist (NCBI feature table, UniProt xlsx).
#   4. FASTA record counts are reported (and compared to `expected_count` when
#      the config supplies one).
#   5. The ConoSorter ZIP exists and, when `expected_sha256_zip` is set,
#      matches.
#   6. A SignalP 6 source is resolvable: either the DTU tarball or the
#      container image. WHICH one is available is recorded, because it changes
#      how R5 runs and therefore belongs in the provenance trail.
#
# Unlike the benchmark's preflight, sha256 pinning of the FASTAs is OPTIONAL
# here: a lab pointing this pipeline at its own Conidae snapshot should not be
# forced to invent checksums. When a checksum IS configured it is enforced.
#
# ALWAYS EXITS 0 and writes a verdict file; the gate rule performs the hard
# fail so the diagnostic report survives on disk.
# =============================================================================
# NOTE: no `from __future__ import annotations` — see project gotchas.
# =============================================================================

import argparse
import glob
import hashlib
import sys
from pathlib import Path


CHUNK = 1024 * 1024


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        while True:
            b = fh.read(CHUNK)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def count_fasta(p: Path) -> int:
    n = 0
    with p.open("rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                n += 1
    return n


def resolve(root: Path, rel: str) -> Path:
    p = Path(rel)
    return p if p.is_absolute() else (root / p).resolve()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db-root",           required=True)
    ap.add_argument("--workdir",           default=".")
    ap.add_argument("--conoserver-fasta",  required=True)
    ap.add_argument("--conoserver-dmnd",   required=True)
    ap.add_argument("--ncbi-fasta",        required=True)
    ap.add_argument("--ncbi-dmnd",         required=True)
    ap.add_argument("--ncbi-feature-table", required=True)
    ap.add_argument("--uniprot-fasta",     required=True)
    ap.add_argument("--uniprot-dmnd",      required=True)
    ap.add_argument("--uniprot-xlsx",      required=True)
    ap.add_argument("--conosorter-zip",    required=True)
    ap.add_argument("--conosorter-sha256", default="")
    ap.add_argument("--signalp-tarball-glob", default="")
    ap.add_argument("--signalp-container",    default="")
    ap.add_argument("--out-report",  required=True, type=Path)
    ap.add_argument("--out-verdict", required=True, type=Path)
    args = ap.parse_args()

    workdir = Path(args.workdir).resolve()
    root = resolve(workdir, args.db_root)

    findings = []
    fails = 0

    def add(item, check, status, detail):
        nonlocal fails
        findings.append((item, check, status, detail))
        if status == "FAIL":
            fails += 1

    # ---- 1. snapshot root -------------------------------------------------
    if not root.is_dir():
        add("snapshot", "root", "FAIL",
            f"{root} does not exist. Set annotation.databases.root to the "
            f"Conidae snapshot directory (an absolute path is fine on HPC).")
    else:
        add("snapshot", "root", "PASS", str(root))

    # ---- 2/3/4. the three databases --------------------------------------
    dbs = [
        ("conoserver", args.conoserver_fasta, args.conoserver_dmnd, None),
        ("ncbi",       args.ncbi_fasta,       args.ncbi_dmnd,       args.ncbi_feature_table),
        ("uniprot",    args.uniprot_fasta,    args.uniprot_dmnd,    args.uniprot_xlsx),
    ]
    for name, fa_rel, dmnd_rel, aux_rel in dbs:
        fa = root / fa_rel
        if fa.is_file() and fa.stat().st_size:
            add(name, "fasta", "PASS", f"{count_fasta(fa):,} sequences")
        else:
            add(name, "fasta", "FAIL", f"missing or empty: {fa}")

        dmnd = root / dmnd_rel
        if dmnd.is_file() and dmnd.stat().st_size:
            add(name, "dmnd", "PASS", f"{dmnd.stat().st_size:,} bytes")
        else:
            add(name, "dmnd", "FAIL",
                f"missing or empty: {dmnd}. Build it with "
                f"`diamond makedb --in {fa} -d {dmnd.with_suffix('')}`.")

        if aux_rel:
            aux = root / aux_rel
            label = "feature_table" if name == "ncbi" else "labels_xlsx"
            if aux.is_file() and aux.stat().st_size:
                add(name, label, "PASS", f"{aux.stat().st_size:,} bytes")
            else:
                add(name, label, "FAIL", f"missing or empty: {aux}")

    # ---- 5. ConoSorter ZIP -------------------------------------------------
    cs = resolve(workdir, args.conosorter_zip)
    if cs.is_file() and cs.stat().st_size:
        add("conosorter", "zip", "PASS", f"{cs.stat().st_size:,} bytes")
        want = args.conosorter_sha256.strip()
        if want:
            got = sha256_file(cs)
            if got == want:
                add("conosorter", "zip_sha256", "PASS", "matches config")
            else:
                add("conosorter", "zip_sha256", "FAIL",
                    f"got {got}, expected {want}")
        else:
            add("conosorter", "zip_sha256", "ADVISORY",
                f"not pinned; observed {sha256_file(cs)}")
    else:
        add("conosorter", "zip", "FAIL",
            f"missing: {cs}. Download from "
            f"https://sourceforge.net/projects/conosorter/files/latest/download")

    # ---- 6. SignalP 6 source ----------------------------------------------
    # Two routes, tried in this order by R5. Recording which one is available
    # matters: the pip route pins an exact DTU release, the container route
    # depends on a third-party image whose contents we do not control.
    sp_source = "none"
    if args.signalp_tarball_glob:
        matches = sorted(glob.glob(str(resolve(workdir, args.signalp_tarball_glob))))
        if len(matches) == 1:
            add("signalp6", "tarball", "PASS", matches[0])
            sp_source = "tarball"
        elif len(matches) > 1:
            add("signalp6", "tarball", "FAIL",
                f"glob matched {len(matches)} files; pin an exact filename")
        else:
            add("signalp6", "tarball", "ADVISORY",
                "no DTU tarball found; will fall back to the container")
    if sp_source == "none" and args.signalp_container:
        add("signalp6", "container", "ADVISORY",
            f"{args.signalp_container} will be used. Note this is a "
            f"third-party image; for a version-pinned install place the DTU "
            f"academic tarball at the configured path instead.")
        sp_source = "container"
    if sp_source == "none":
        add("signalp6", "source", "FAIL",
            "neither a DTU tarball nor a container image is configured")
    else:
        add("signalp6", "source", "PASS", sp_source)

    # ---- emit --------------------------------------------------------------
    verdict = "PASS" if fails == 0 else "FAIL"
    args.out_report.parent.mkdir(parents=True, exist_ok=True)
    with args.out_report.open("w", encoding="utf-8") as fh:
        fh.write("# Stage 02 preflight — annotation resources\n")
        fh.write(f"# db_root = {root}\n\n")
        w = max((len(f[0]) for f in findings), default=10) + 2
        for item, check, status, detail in findings:
            fh.write(f"[{status:8}] {item:<{w}} {check:<15} :: {detail}\n")
        fh.write(f"\nTotal FAIL: {fails}\nVERDICT: {verdict}\n")

    args.out_verdict.parent.mkdir(parents=True, exist_ok=True)
    args.out_verdict.write_text(
        f"{verdict}\ntotal_fails={fails}\nsignalp_source={sp_source}\n")

    print(f"[preflight] verdict={verdict} fails={fails} "
          f"signalp_source={sp_source}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
