# =============================================================================
# 00_preprocessing.smk  -  Stage 00: read QC and trimming with fastp
# =============================================================================
# CONTRACT
#   IN   validated sample sheet (sample_id, fastq_1[, fastq_2], layout)
#   OUT  results/{run_id}/00_preprocessing/{sample}/
#          {sample}_clean_R1.fastq.gz  [+ _clean_R2.fastq.gz when PE]
#          {sample}_fastp.html / .json
#          {sample}_fastp.summary.tsv
#        results/{run_id}/00_preprocessing/preprocessing_manifest.tsv (+ .sha256)
#
# TOGGLE. `preprocessing.fastp: false` skips the fastp rule entirely; the
# manifest is still produced, recording `reads_source = raw` and pointing at
# the untouched input FASTQs. Stage 01 reads its inputs through the
# `reads_for()` helper in the Snakefile, so the toggle re-points the whole DAG
# without any rule in Stage 01 knowing about it.
#
# SE/PE. One rule serves both layouts. Snakemake cannot express "this output
# exists only for some wildcard values", so the R2 output is declared for
# every sample and the wrapper writes an empty placeholder for SE libraries.
# The placeholder is never consumed: `reads_for()` returns a single-element
# dict for SE samples, so no assembler rule ever names it.
# =============================================================================


# ---- fastp flag assembly ----------------------------------------------------
def _fastp_extra_flags(sample_id: str):
    """Translate the config block into fastp CLI flags.

    Kept as a function rather than a static string because a few flags are
    layout-dependent (`--detect_adapter_for_pe` is meaningless for SE and
    fastp warns about it) and because adapter overrides are optional.
    """
    c = CFG_PREP
    flags = [
        "--qualified_quality_phred", str(c["qualified_quality_phred"]),
        "--unqualified_percent_limit", str(c["unqualified_percent_limit"]),
        "--n_base_limit", str(c["n_base_limit"]),
        "--length_required", str(c["length_required"]),
    ]
    if c.get("cut_front", True):
        flags.append("--cut_front")
    if c.get("cut_tail", True):
        flags.append("--cut_tail")
    flags += [
        "--cut_window_size", str(c.get("cut_window_size", 4)),
        "--cut_mean_quality", str(c.get("cut_mean_quality", 20)),
    ]
    if is_paired(sample_id) and c.get("detect_adapter_for_pe", True):
        flags.append("--detect_adapter_for_pe")
    if c.get("adapter_sequence"):
        flags += ["--adapter_sequence", c["adapter_sequence"]]
    if is_paired(sample_id) and c.get("adapter_sequence_r2"):
        flags += ["--adapter_sequence_r2", c["adapter_sequence_r2"]]
    if c.get("dedup", False):
        flags.append("--dedup")
    if c.get("extra"):
        flags += str(c["extra"]).split()
    return flags


def _fastp_inputs(wildcards):
    """Raw FASTQs for this sample, straight from the sheet."""
    rec = sample_record(wildcards.sample)
    out = {"r1": rec["fastq_1"]}
    if rec["layout"] == "PE":
        out["r2"] = rec["fastq_2"]
    return out


# ============================================================================
# Stage 00 R1 — fastp
# ============================================================================
rule fastp:
    input:
        unpack(_fastp_inputs),
        sheet = SAMPLES_VALIDATED,
    output:
        r1      = f"{D_PREP}/{{sample}}/{{sample}}_clean_R1.fastq.gz",
        r2      = f"{D_PREP}/{{sample}}/{{sample}}_clean_R2.fastq.gz",
        html    = f"{D_PREP}/{{sample}}/{{sample}}_fastp.html",
        json    = f"{D_PREP}/{{sample}}/{{sample}}_fastp.json",
        summary = f"{D_PREP}/{{sample}}/{{sample}}_fastp.summary.tsv",
    params:
        layout = lambda w: "PE" if is_paired(w.sample) else "SE",
        in_r2  = lambda w: (sample_record(w.sample)["fastq_2"]
                            if is_paired(w.sample) else ""),
        # For SE the wrapper is handed an empty out_r2 and skips writing it;
        # the rule then touches the declared placeholder (see shell below).
        out_r2 = lambda w, output: output.r2 if is_paired(w.sample) else "",
        extra  = lambda w: " ".join(
            f'"{f}"' for f in _fastp_extra_flags(w.sample)),
    threads:
        CFG_PREP.get("threads", 4)
    resources:
        mem_mb = CFG_RES.get("fastp_mem_mb", 8192),
    conda:
        "../envs/fastp.yaml"
    log:
        "logs/00_preprocessing/fastp__{sample}.log"
    shell:
        r"""
        bash workflow/scripts/01_run_fastp.sh \
            "{params.layout}" \
            "{input.r1}" \
            "{params.in_r2}" \
            "{output.r1}" \
            "{params.out_r2}" \
            "{output.html}" \
            "{output.json}" \
            "{output.summary}" \
            "{threads}" \
            "{wildcards.sample}" \
            {params.extra} \
            > {log} 2>&1

        # Single-end libraries have no R2. Snakemake requires every declared
        # output to exist, so materialise an empty placeholder. It is never
        # consumed: reads_for() returns only r1 for SE samples.
        if [ "{params.layout}" = "SE" ]; then
            : | gzip -c > "{output.r2}"
        fi
        """


# ============================================================================
# Stage 00 R2 — concatenate a group's libraries into one read set
# ----------------------------------------------------------------------------
# The assembly unit is the GROUP, so the libraries that share a `group` value
# are pooled here, AFTER each has been cleaned individually. Cleaning first and
# pooling second is deliberate: adapter content and quality vary between lanes
# and sequencing runs, so per-library detection is more accurate, and a bad
# replicate stays visible in its own QC report instead of being averaged away.
#
# Concatenation order is the sorted sample_id order fixed by
# `libraries_in_group()`. Assemblers are order-sensitive, so pinning it keeps
# runs reproducible.
#
# Groups of ONE library still pass through this rule. The extra copy costs
# disk but buys a uniform input path for Stage 01, which then needs no
# single-vs-pooled branching anywhere.
#
# gzip members concatenate losslessly at the byte level -- `cat a.gz b.gz` is a
# valid gzip stream -- so no decompression round-trip is needed.
# ============================================================================
def _group_concat_inputs(wildcards):
    members = libraries_in_group(wildcards.group)
    r1 = [reads_for(m["sample_id"])["r1"] for m in members]
    ins = {"r1": r1}
    if group_is_paired(wildcards.group):
        ins["r2"] = [reads_for(m["sample_id"])["r2"] for m in members]
    return ins


rule concat_group_reads:
    input:
        unpack(_group_concat_inputs),
    output:
        r1      = f"{D_PREP}/_groups/{{group}}/{{group}}_R1.fastq.gz",
        r2      = f"{D_PREP}/_groups/{{group}}/{{group}}_R2.fastq.gz",
        summary = f"{D_PREP}/_groups/{{group}}/{{group}}_concat.summary.tsv",
    params:
        layout    = lambda w: group_layout(w.group),
        n_libs    = lambda w: len(libraries_in_group(w.group)),
        lib_ids   = lambda w: ",".join(m["sample_id"]
                                       for m in libraries_in_group(w.group)),
        r2_inputs = lambda w, input: (" ".join(input.r2)
                                      if group_is_paired(w.group) else ""),
    resources:
        mem_mb = CFG_RES.get("py_mem_mb", 4096),
    conda:
        "../envs/fastp.yaml"          # provides seqkit for the read counts
    log:
        "logs/00_preprocessing/concat__{group}.log"
    shell:
        r"""
        exec > >(tee -a "{log}") 2>&1
        set -euo pipefail

        mkdir -p "$(dirname "{output.r1}")"

        echo "[concat] group={wildcards.group} layout={params.layout} \
libraries={params.n_libs} ({params.lib_ids})"

        # Concatenating gzip members is byte-safe: the result is a valid gzip
        # stream, so no decompress/recompress cycle is required.
        cat {input.r1} > "{output.r1}"

        if [ "{params.layout}" = "PE" ]; then
            cat {params.r2_inputs} > "{output.r2}"
        else
            # Declared output must exist; never consumed for SE groups because
            # reads_for_group() returns only r1.
            : | gzip -c > "{output.r2}"
        fi

        # Verify the concatenation produced readable gzip streams. A truncated
        # member here would silently shrink every downstream assembly.
        gzip -t "{output.r1}"
        [ "{params.layout}" = "PE" ] && gzip -t "{output.r2}"

        N_R1=$(seqkit stats -T "{output.r1}" 2>/dev/null | awk 'NR==2{{print $4}}')
        N_R2="NA"
        if [ "{params.layout}" = "PE" ]; then
            N_R2=$(seqkit stats -T "{output.r2}" 2>/dev/null | awk 'NR==2{{print $4}}')
            if [ "$N_R1" != "$N_R2" ]; then
                echo "[concat] ERROR: R1/R2 record counts differ ($N_R1 vs $N_R2);" \
                     "the pooled libraries are not consistently paired." >&2
                exit 3
            fi
        fi

        {{
            printf 'metric\tvalue\n'
            printf 'group\t%s\n'          "{wildcards.group}"
            printf 'layout\t%s\n'         "{params.layout}"
            printf 'n_libraries\t%s\n'    "{params.n_libs}"
            printf 'library_ids\t%s\n'    "{params.lib_ids}"
            printf 'reads_r1\t%s\n'       "$N_R1"
            printf 'reads_r2\t%s\n'       "$N_R2"
            printf 'concat_r1\t%s\n'      "{output.r1}"
            printf 'concat_r2\t%s\n'      "{output.r2}"
        }} > "{output.summary}"

        echo "[concat] done: $N_R1 records in R1"
        """


# ============================================================================
# Stage 00 R3 — manifest
# ----------------------------------------------------------------------------
# Runs in both modes. When fastp is enabled it depends on every per-sample
# summary; when disabled it depends only on the validated sheet, so the rule
# is cheap and the results tree keeps the same shape either way.
# ============================================================================
def _manifest_inputs(wildcards=None):
    ins = {"samples": SAMPLES_VALIDATED}
    if CFG_PREP.get("fastp", True):
        ins["summaries"] = [
            f"{D_PREP}/{s}/{s}_fastp.summary.tsv" for s in sample_ids()
        ]
    # Pooling summaries: depending on them here makes the manifest the single
    # completion marker for the whole of Stage 00, so `rule preprocessing`
    # cannot report success while a group concatenation is still missing.
    ins["concats"] = [
        f"{D_PREP}/_groups/{g}/{g}_concat.summary.tsv" for g in group_ids()
    ]
    return ins


rule preprocessing_manifest:
    input:
        unpack(_manifest_inputs),
    output:
        manifest = f"{D_PREP}/preprocessing_manifest.tsv",
        sha256   = f"{D_PREP}/preprocessing_manifest.tsv.sha256",
    params:
        summary_dir   = D_PREP,
        group_dir     = f"{D_PREP}/_groups",
        fastp_enabled = "true" if CFG_PREP.get("fastp", True) else "false",
    resources:
        mem_mb = CFG_RES.get("py_mem_mb", 4096),
    conda:
        "../envs/python_base.yaml"
    log:
        "logs/00_preprocessing/preprocessing_manifest.log"
    shell:
        # `python`, not `python3`: inside a conda env the former is guaranteed
        # to be the env interpreter.
        r"""
        python workflow/scripts/02_build_preprocessing_manifest.py \
            --samples       "{input.samples}" \
            --summary-dir   "{params.summary_dir}" \
            --group-dir     "{params.group_dir}" \
            --fastp-enabled "{params.fastp_enabled}" \
            --out-manifest  "{output.manifest}" \
            > {log} 2>&1
        """


# ---- Stage targets ----------------------------------------------------------
def stage00_targets():
    """Files activated by Stage 00 (consumed by `rule all` and `rule
    preprocessing`)."""
    return [
        f"{D_PREP}/preprocessing_manifest.tsv",
        f"{D_PREP}/preprocessing_manifest.tsv.sha256",
    ]


localrules: preprocessing_manifest
