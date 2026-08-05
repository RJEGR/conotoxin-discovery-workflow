# =============================================================================
# 01_assembly.smk  -  Stage 01: de novo assembly + reference-free filtering
# =============================================================================
# THE UNIT OF THIS STAGE IS THE GROUP, NOT THE LIBRARY.
#
# Stage 00 cleans each library separately and then concatenates the libraries
# that share a `group` value. Stage 01 assembles that pooled read set once.
# Assembling replicates independently fragments transcripts that the combined
# depth recovers whole — a disproportionate loss for conotoxins, whose
# precursors are short and unevenly expressed, so a partial ORF from one
# replicate typically fails the Stage 02 signal-peptide test and disappears
# from the annotation.
#
# A library with no `group` value forms a group of one, which reproduces
# one-assembly-per-library behaviour.
#
# CONTRACT
#   IN   reads_for_group(group) -> the concatenated read set from Stage 00
#   OUT  results/{run_id}/01_assembly/{group}/{asm}/
#          contigs.fasta            raw assembly
#          run.summary.tsv          assembler metrics
#          transrate/               TransRate outputs (paired-end only)
#          transrate.summary.tsv    filtering metrics
#          filtered_contigs.fasta   <- THE input Stage 02 consumes
#        results/{run_id}/01_assembly/assembly_manifest.tsv (+ .sha256)
#
# THE filtered_contigs.fasta CONTRACT. Stage 02 never branches on layout: it
# always reads filtered_contigs.fasta. That file is TransRate's good.*.fa for
# paired-end groups and a copy of contigs.fasta for single-end ones, since
# TransRate's read-metrics mode requires both mates. `filter_method` in the
# manifest records which happened; SE and PE contig counts are not comparable
# without stratifying on it.
#
# CONTAINERISATION. Trinity and TransRate run in containers via
# container_run.sh (Docker on macOS, Apptainer on SLURM, auto-detected).
# rnaSPAdes and RNA-Bloom run natively from conda: SPAdes has osx-64 builds
# that work under Rosetta, RNA-Bloom is noarch Java. Fewer containers means
# fewer moving parts on HPC.
# =============================================================================


# ---- Shared helpers ---------------------------------------------------------
def _asm_cfg(asm: str):
    return CFG_ASM["assemblers"][asm]


def _assembly_reads(wildcards):
    """The group's pooled reads, produced by Stage 00's concat rule."""
    return reads_for_group(wildcards.group)


def _r2_or_empty(wildcards, input_obj) -> str:
    """R2 path when the group is paired, empty string otherwise.

    Read off `input` rather than recomputed, so the wrapper argument and the
    declared dependency cannot drift apart.
    """
    return (getattr(input_obj, "r2", "")
            if group_is_paired(wildcards.group) else "")


# Container arguments, shared by the Trinity and TransRate rules.
_CT_RUNTIME  = CFG_CONTAINERS["runtime"]
_CT_PLATFORM = CFG_CONTAINERS.get("docker_platform", "")
_CT_CACHE    = CFG_CONTAINERS.get("cache_dir", "resources/containers")
_CT_BINDS    = ",".join(CFG_CONTAINERS.get("extra_binds", []))


# ============================================================================
# Stage 01 R1a — Trinity  (containerised)
# ============================================================================
rule assemble_trinity:
    input:
        unpack(_assembly_reads),
    output:
        contigs = f"{D_ASM}/{{group}}/trinity/contigs.fasta",
        summary = f"{D_ASM}/{{group}}/trinity/run.summary.tsv",
    params:
        layout       = lambda w: group_layout(w.group),
        r2           = lambda w, input: _r2_or_empty(w, input),
        mem_gb       = _asm_cfg("trinity")["max_memory_gb"],
        no_normalize = "true" if _asm_cfg("trinity").get("no_normalize_reads", False) else "false",
        image        = _asm_cfg("trinity")["container"],
        runtime      = _CT_RUNTIME,
        platform     = _CT_PLATFORM,
        cache_dir    = _CT_CACHE,
        binds        = _CT_BINDS,
    threads:
        _asm_cfg("trinity")["threads"]
    resources:
        mem_mb = CFG_RES.get("assembly_mem_mb", 32768),
    conda:
        "../envs/container_host.yaml"
    log:
        "logs/01_assembly/trinity__{group}.log"
    shell:
        r"""
        bash workflow/scripts/03_run_trinity.sh \
            "{params.layout}" "{input.r1}" "{params.r2}" \
            "{output.contigs}" "{output.summary}" \
            "{params.mem_gb}" "{threads}" "{wildcards.group}" \
            "{params.no_normalize}" \
            "{params.runtime}" "{params.image}" "{params.platform}" \
            "{params.cache_dir}" "{params.binds}" \
            > {log} 2>&1
        """


# ============================================================================
# Stage 01 R1b — rnaSPAdes  (native conda)
# ============================================================================
rule assemble_rnaspades:
    input:
        unpack(_assembly_reads),
    output:
        contigs = f"{D_ASM}/{{group}}/rnaspades/contigs.fasta",
        summary = f"{D_ASM}/{{group}}/rnaspades/run.summary.tsv",
    params:
        layout = lambda w: group_layout(w.group),
        r2     = lambda w, input: _r2_or_empty(w, input),
        mem_gb = _asm_cfg("rnaspades")["max_memory_gb"],
    threads:
        _asm_cfg("rnaspades")["threads"]
    resources:
        mem_mb = CFG_RES.get("assembly_mem_mb", 32768),
    conda:
        "../envs/assembly_rnaspades.yaml"
    log:
        "logs/01_assembly/rnaspades__{group}.log"
    shell:
        r"""
        bash workflow/scripts/04_run_rnaspades.sh \
            "{params.layout}" "{input.r1}" "{params.r2}" \
            "{output.contigs}" "{output.summary}" \
            "{params.mem_gb}" "{threads}" "{wildcards.group}" \
            > {log} 2>&1
        """


# ============================================================================
# Stage 01 R1c — RNA-Bloom  (native conda, noarch Java)
# ============================================================================
rule assemble_rnabloom:
    input:
        unpack(_assembly_reads),
    output:
        contigs = f"{D_ASM}/{{group}}/rnabloom/contigs.fasta",
        summary = f"{D_ASM}/{{group}}/rnabloom/run.summary.tsv",
    params:
        layout = lambda w: group_layout(w.group),
        r2     = lambda w, input: _r2_or_empty(w, input),
        mem_gb = _asm_cfg("rnabloom")["max_memory_gb"],
    threads:
        _asm_cfg("rnabloom")["threads"]
    resources:
        mem_mb = CFG_RES.get("assembly_mem_mb", 32768),
    conda:
        "../envs/assembly_rnabloom.yaml"
    log:
        "logs/01_assembly/rnabloom__{group}.log"
    shell:
        r"""
        bash workflow/scripts/05_run_rnabloom.sh \
            "{params.layout}" "{input.r1}" "{params.r2}" \
            "{output.contigs}" "{output.summary}" \
            "{params.mem_gb}" "{threads}" "{wildcards.group}" \
            > {log} 2>&1
        """


# ============================================================================
# Stage 01 R2 — TransRate filtering / filtered_contigs.fasta
# ----------------------------------------------------------------------------
# One rule for all three assemblers (the `asm` wildcard). Emits
# filtered_contigs.fasta unconditionally, so the Stage 02 input path is the
# same for every group and every assembler regardless of layout.
# ============================================================================
def _transrate_inputs(wildcards):
    ins = dict(reads_for_group(wildcards.group))
    ins["contigs"] = (f"{D_ASM}/{wildcards.group}/{wildcards.asm}/"
                      f"contigs.fasta")
    return ins


rule transrate_filter:
    input:
        unpack(_transrate_inputs),
    output:
        filtered = f"{D_ASM}/{{group}}/{{asm}}/filtered_contigs.fasta",
        summary  = f"{D_ASM}/{{group}}/{{asm}}/transrate.summary.tsv",
    params:
        layout    = lambda w: group_layout(w.group),
        r2        = lambda w, input: _r2_or_empty(w, input),
        tr_dir    = lambda w: f"{D_ASM}/{w.group}/{w.asm}/transrate",
        enabled   = "true" if CFG_ASM["transrate"].get("enabled", True) else "false",
        sanitize  = "true" if CFG_ASM["transrate"].get("sanitize_reads", True) else "false",
        image     = CFG_ASM["transrate"]["container"],
        runtime   = _CT_RUNTIME,
        platform  = _CT_PLATFORM,
        cache_dir = _CT_CACHE,
        binds     = _CT_BINDS,
    threads:
        CFG_ASM["transrate"].get("threads", 4)
    resources:
        mem_mb = CFG_RES.get("transrate_mem_mb", 8192),
    conda:
        "../envs/container_host.yaml"
    log:
        "logs/01_assembly/transrate__{group}__{asm}.log"
    shell:
        r"""
        bash workflow/scripts/06_run_transrate.sh \
            "{params.layout}" "{input.contigs}" "{input.r1}" "{params.r2}" \
            "{output.filtered}" "{output.summary}" "{params.tr_dir}" \
            "{threads}" "{wildcards.group}" "{wildcards.asm}" \
            "{params.enabled}" "{params.sanitize}" \
            "{params.runtime}" "{params.image}" "{params.platform}" \
            "{params.cache_dir}" "{params.binds}" \
            > {log} 2>&1
        """


# ============================================================================
# Stage 01 R3 — manifest
# ============================================================================
def _assembly_manifest_inputs(wildcards=None):
    outs = []
    for g in group_ids():
        for asm in enabled_assemblers():
            outs.append(f"{D_ASM}/{g}/{asm}/filtered_contigs.fasta")
            outs.append(f"{D_ASM}/{g}/{asm}/transrate.summary.tsv")
            outs.append(f"{D_ASM}/{g}/{asm}/run.summary.tsv")
    return outs


rule assembly_manifest:
    input:
        groups = GROUPS_VALIDATED,
        parts  = _assembly_manifest_inputs,
    output:
        manifest = f"{D_ASM}/assembly_manifest.tsv",
        sha256   = f"{D_ASM}/assembly_manifest.tsv.sha256",
    params:
        assembly_dir = D_ASM,
        assemblers   = ",".join(enabled_assemblers()),
        min_contigs  = CFG_ASM["invariants"].get("min_contigs_per_assembly", 1),
    resources:
        mem_mb = CFG_RES.get("py_mem_mb", 4096),
    conda:
        "../envs/python_base.yaml"
    log:
        "logs/01_assembly/assembly_manifest.log"
    shell:
        r"""
        python workflow/scripts/07_build_assembly_manifest.py \
            --groups       "{input.groups}" \
            --assembly-dir "{params.assembly_dir}" \
            --assemblers   "{params.assemblers}" \
            --min-contigs  "{params.min_contigs}" \
            --out-manifest "{output.manifest}" \
            > {log} 2>&1
        """


# ---- Stage targets ----------------------------------------------------------
def stage01_targets():
    """Files activated by Stage 01 (consumed by `rule all` / `rule assembly`)."""
    return [
        f"{D_ASM}/assembly_manifest.tsv",
        f"{D_ASM}/assembly_manifest.tsv.sha256",
    ]


localrules: assembly_manifest
