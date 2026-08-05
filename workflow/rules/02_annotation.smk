# =============================================================================
# 02_annotation.smk  -  Stage 02: conotoxin annotation and confidence tiers
# =============================================================================
# THE UNIT OF THIS STAGE IS THE GROUP. Stage 01 assembled each group's pooled
# reads with every enabled assembler; here those assemblies are dereplicated
# into one non-redundant peptide set per group and annotated once.
#
# DAG (per group unless noted)
#   R1  preflight + gate            (1x)  databases, ConoSorter ZIP, SignalP source
#   R2  predict_orfs          (per asm)   TransDecoder LongOrfs + DIAMOND + Predict
#   R3  dedup_merge_pep                   MMseqs2 easy-cluster over the assemblers
#   R4a conosorter_run                    ConoSorter RegEx + pHMM in one call
#   R4b conosorter_merge                  outer join -> merge.tab / merge.fasta
#   R5  signalp_and_prepare               SignalP 6 + filter + split
#   R6a diamond_annotate                  blastp x 3 DBs x 2 parts
#   R6b parse_besthits                    outfmt6 -> best hits + coverage audit
#   R7  annotate_labels                   ConoServer / NCBI / UniProt labels
#   R8  normalize_final                   vendored legacy-11 normalisation
#   R9  classify_confidence               vendored legacy-12 decision tree
#   R10 aggregate + invariants + gate (1x)
#
# `no_refine_starts` is FALSE here (it is true in the simulation benchmark).
# Real libraries assemble far more contigs than TransDecoder's start-codon PWM
# needs to train, and refinement sharpens exactly the region this pipeline is
# most sensitive to: the ORF N-terminus IS the signal peptide that R5 feeds to
# SignalP, and precursors are short enough that a few misplaced codons shift a
# large fraction of the molecule.
# =============================================================================

# ---- Config handles ---------------------------------------------------------
ANN_DB   = CFG_ANN["databases"]
ANN_CS   = CFG_ANN["conosorter"]
ANN_SP   = CFG_ANN["signalp6"]
ANN_PRED = CFG_ANN["predict"]
ANN_DED  = CFG_ANN["dedup"]
ANN_DIA  = CFG_ANN["diamond_annotate"]
ANN_BH   = CFG_ANN["besthits"]
ANN_INV  = CFG_ANN["invariants"]

DB_ROOT = ANN_DB["root"]


def _db(key, field):
    return f"{DB_ROOT}/{ANN_DB[key][field]}"


# ============================================================================
# R1 — preflight + gate
# ============================================================================
rule annotation_preflight:
    output:
        report  = f"{D_QC}/annotation_preflight.txt",
        verdict = f"{D_QC}/annotation_preflight_verdict.txt",
    params:
        db_root   = DB_ROOT,
        cs_fa     = ANN_DB["conoserver"]["fasta"],
        cs_dmnd   = ANN_DB["conoserver"]["dmnd"],
        nc_fa     = ANN_DB["ncbi"]["fasta"],
        nc_dmnd   = ANN_DB["ncbi"]["dmnd"],
        nc_ft     = ANN_DB["ncbi"]["feature_table"],
        up_fa     = ANN_DB["uniprot"]["fasta"],
        up_dmnd   = ANN_DB["uniprot"]["dmnd"],
        up_xlsx   = ANN_DB["uniprot"]["labels_xlsx"],
        cs_zip    = ANN_CS["local_zip"],
        cs_sha    = ANN_CS.get("expected_sha256_zip", ""),
        sp_tar    = ANN_SP.get("tarball_glob", ""),
        sp_cont   = ANN_SP.get("container", ""),
    conda:
        "../envs/python_base.yaml"
    log:
        "logs/02_annotation/preflight.log"
    shell:
        r"""
        python workflow/scripts/09_preflight_dbs.py \
            --db-root "{params.db_root}" --workdir "." \
            --conoserver-fasta "{params.cs_fa}" --conoserver-dmnd "{params.cs_dmnd}" \
            --ncbi-fasta "{params.nc_fa}" --ncbi-dmnd "{params.nc_dmnd}" \
            --ncbi-feature-table "{params.nc_ft}" \
            --uniprot-fasta "{params.up_fa}" --uniprot-dmnd "{params.up_dmnd}" \
            --uniprot-xlsx "{params.up_xlsx}" \
            --conosorter-zip "{params.cs_zip}" --conosorter-sha256 "{params.cs_sha}" \
            --signalp-tarball-glob "{params.sp_tar}" \
            --signalp-container "{params.sp_cont}" \
            --out-report "{output.report}" --out-verdict "{output.verdict}" \
            > {log} 2>&1
        """


rule annotation_preflight_gate:
    input:
        verdict = f"{D_QC}/annotation_preflight_verdict.txt",
        report  = f"{D_QC}/annotation_preflight.txt",
    output:
        ok = f"{D_QC}/annotation_preflight.OK",
    shell:
        # Both files are declared as inputs so the report path can be quoted
        # directly. Bash parameter expansion (${var%suffix}) cannot be used
        # here: Snakemake parses `{...}` as its own format field first and the
        # rule fails at parse time.
        r"""
        v=$(head -1 {input.verdict})
        if [ "$v" != "PASS" ]; then
            echo "Stage 02 preflight verdict=$v — halting." >&2
            echo "See {input.report} for the failing checks." >&2
            exit 1
        fi
        touch {output.ok}
        """


# ============================================================================
# R2 — ORF prediction (one job per group x assembler)
# ============================================================================
rule predict_orfs:
    input:
        contigs = f"{D_ASM}/{{group}}/{{asm}}/filtered_contigs.fasta",
        db_dmnd = _db("uniprot", "dmnd"),
        gate    = f"{D_QC}/annotation_preflight.OK",
    output:
        pep     = f"{D_ANN}/{{group}}/{{asm}}/predict/transdecoder.pep",
        cds     = f"{D_ANN}/{{group}}/{{asm}}/predict/transdecoder.cds",
        gff3    = f"{D_ANN}/{{group}}/{{asm}}/predict/transdecoder.gff3",
        bed     = f"{D_ANN}/{{group}}/{{asm}}/predict/transdecoder.bed",
        blast   = f"{D_ANN}/{{group}}/{{asm}}/predict/blastp.support.outfmt6",
        summary = f"{D_ANN}/{{group}}/{{asm}}/predict/predict.summary.tsv",
    params:
        outdir   = lambda w: f"{D_ANN}/{w.group}/{w.asm}/predict",
        min_orf  = ANN_PRED["min_orf_length"],
        evalue   = ANN_PRED["support_evalue"],
        mts      = ANN_PRED["support_max_target_seqs"],
        single   = "true" if ANN_PRED.get("single_best_only", False) else "false",
        norefine = "true" if ANN_PRED.get("no_refine_starts", False) else "false",
    threads:
        ANN_PRED.get("threads", 8)
    resources:
        mem_mb = CFG_RES.get("predict_mem_mb", 16384),
    conda:
        "../envs/annotation.yaml"
    log:
        "logs/02_annotation/predict__{group}__{asm}.log"
    shell:
        r"""
        bash workflow/scripts/10_predict_orfs.sh \
            "{input.contigs}" "{input.db_dmnd}" "{params.outdir}" \
            "{params.min_orf}" "{params.evalue}" "{params.mts}" \
            "{threads}" "{params.single}" "{params.norefine}" \
            > {log} 2>&1
        """


# ============================================================================
# R3 — cross-assembler dereplication
# ============================================================================
def _peps_for_group(wildcards):
    return [f"{D_ANN}/{wildcards.group}/{asm}/predict/transdecoder.pep"
            for asm in enabled_assemblers()]


rule dedup_merge_pep:
    input:
        peps = _peps_for_group,
    output:
        nr_pep   = f"{D_ANN}/{{group}}/dedup/nr.pep.fasta",
        clusters = f"{D_ANN}/{{group}}/dedup/nr.cluster.tsv",
        combined = f"{D_ANN}/{{group}}/dedup/combined.pep.fasta",
        audit    = f"{D_ANN}/{{group}}/dedup/dedup.audit.tsv",
        summary  = f"{D_ANN}/{{group}}/dedup/dedup.summary.tsv",
    params:
        outdir = lambda w: f"{D_ANN}/{w.group}/dedup",
        minid  = ANN_DED["min_seq_id"],
        cov    = ANN_DED["coverage"],
        covmod = ANN_DED["cov_mode"],
        sens   = ANN_DED["sensitivity"],
    threads:
        ANN_DED.get("threads", 8)
    resources:
        mem_mb = CFG_RES.get("dedup_mem_mb", 16384),
    conda:
        "../envs/annotation.yaml"
    log:
        "logs/02_annotation/dedup__{group}.log"
    shell:
        r"""
        bash workflow/scripts/11_dedup_merge_pep.sh \
            "{params.outdir}" "{params.minid}" "{params.cov}" \
            "{params.covmod}" "{params.sens}" "{threads}" \
            {input.peps} > {log} 2>&1
        """


# ============================================================================
# R4a/R4b — ConoSorter and merge
# ============================================================================
rule conosorter_run:
    input:
        nr_pep = f"{D_ANN}/{{group}}/dedup/nr.pep.fasta",
        cs_zip = ANN_CS["local_zip"],
    output:
        regex   = f"{D_ANN}/{{group}}/conosorter/Regex.tab",
        phmm    = f"{D_ANN}/{{group}}/conosorter/pHMM.tab",
        summary = f"{D_ANN}/{{group}}/conosorter/run.summary.tsv",
    params:
        outdir = lambda w: f"{D_ANN}/{w.group}/conosorter",
    threads:
        CFG_RES.get("py_threads", 2)
    resources:
        mem_mb = CFG_RES.get("conosorter_mem_mb", 4096),
    conda:
        "../envs/conosorter.yaml"
    log:
        "logs/02_annotation/conosorter__{group}.log"
    shell:
        r"""
        bash workflow/scripts/12_conosorter_run.sh \
            "{input.nr_pep}" "{input.cs_zip}" "{params.outdir}" \
            > {log} 2>&1
        """


rule conosorter_merge:
    input:
        regex = f"{D_ANN}/{{group}}/conosorter/Regex.tab",
        phmm  = f"{D_ANN}/{{group}}/conosorter/pHMM.tab",
    output:
        tab   = f"{D_ANN}/{{group}}/conosorter/merge.tab",
        fasta = f"{D_ANN}/{{group}}/conosorter/merge.fasta",
        audit = f"{D_ANN}/{{group}}/conosorter/merge.audit.tsv",
    resources:
        mem_mb = CFG_RES.get("py_mem_mb", 4096),
    conda:
        "../envs/python_base.yaml"
    log:
        "logs/02_annotation/conosorter_merge__{group}.log"
    shell:
        r"""
        python workflow/scripts/13_merge_conosorter.py \
            --regex-tab "{input.regex}" --phmm-tab "{input.phmm}" \
            --out-tab "{output.tab}" --out-fasta "{output.fasta}" \
            --audit "{output.audit}" > {log} 2>&1
        """


# ============================================================================
# R5 — SignalP 6 + candidate preparation
# ============================================================================
rule signalp_and_prepare:
    input:
        tab   = f"{D_ANN}/{{group}}/conosorter/merge.tab",
        fasta = f"{D_ANN}/{{group}}/conosorter/merge.fasta",
    output:
        pred     = f"{D_ANN}/{{group}}/signalp/signalp_output/prediction_results.txt",
        filtered = f"{D_ANN}/{{group}}/signalp/filtered.tab",
        reorg    = f"{D_ANN}/{{group}}/signalp/reorganized.xlsx",
        pos_xlsx = f"{D_ANN}/{{group}}/signalp/positivo_superfamilia.xlsx",
        neg_xlsx = f"{D_ANN}/{{group}}/signalp/negativo_superfamilia.xlsx",
        pos_fa   = f"{D_ANN}/{{group}}/signalp/positivo.fasta",
        neg_fa   = f"{D_ANN}/{{group}}/signalp/negativo.fasta",
        audit    = f"{D_ANN}/{{group}}/signalp/signalp_prepare.audit.tsv",
    params:
        outdir    = lambda w: f"{D_ANN}/{w.group}/signalp",
        tarball   = ANN_SP.get("tarball_glob", ""),
        container = ANN_SP.get("container", ""),
        organism  = ANN_SP.get("organism", "eukarya"),
        mode      = ANN_SP.get("mode", "fast"),
        fmt       = ANN_SP.get("format", "txt"),
        runtime   = CFG_CONTAINERS["runtime"],
        platform  = CFG_CONTAINERS.get("docker_platform", ""),
        cache     = CFG_CONTAINERS.get("cache_dir", "resources/containers"),
        binds     = ",".join(CFG_CONTAINERS.get("extra_binds", [])),
    threads:
        CFG_RES.get("py_threads", 2)
    resources:
        mem_mb = CFG_RES.get("signalp_mem_mb", 12288),
    conda:
        "../envs/signalp.yaml"
    log:
        "logs/02_annotation/signalp__{group}.log"
    shell:
        r"""
        python workflow/scripts/14_signalp_and_prepare.py \
            --merge-tab "{input.tab}" --merge-fasta "{input.fasta}" \
            --outdir "{params.outdir}" \
            --signalp-tarball "{params.tarball}" \
            --signalp-container "{params.container}" \
            --organism "{params.organism}" --mode "{params.mode}" \
            --format "{params.fmt}" \
            --container-runtime "{params.runtime}" \
            --container-platform "{params.platform}" \
            --container-cache "{params.cache}" \
            --container-binds "{params.binds}" \
            > {log} 2>&1
        """


# ============================================================================
# R6a/R6b — DIAMOND annotation and best-hit selection
# ============================================================================
rule diamond_annotate:
    input:
        pos = f"{D_ANN}/{{group}}/signalp/positivo.fasta",
        neg = f"{D_ANN}/{{group}}/signalp/negativo.fasta",
        db_uniprot    = _db("uniprot", "dmnd"),
        db_ncbi       = _db("ncbi", "dmnd"),
        db_conoserver = _db("conoserver", "dmnd"),
    output:
        outs = expand(
            f"{D_ANN}/{{{{group}}}}/annotate/blast/{{part}}__{{db}}.outfmt6",
            part=["positivo", "negativo"],
            db=["uniprot", "ncbi", "conoserver"]),
        summary = f"{D_ANN}/{{group}}/annotate/blast/diamond.summary.tsv",
    params:
        outdir  = lambda w: f"{D_ANN}/{w.group}/annotate/blast",
        evalue  = ANN_DIA["evalue"],
        qcov    = ANN_DIA["min_qcov"],
        hsps    = ANN_DIA["max_hsps"],
        mts_up  = ANN_DIA["max_target_seqs"]["uniprot"],
        mts_nc  = ANN_DIA["max_target_seqs"]["ncbi"],
        mts_cs  = ANN_DIA["max_target_seqs"]["conoserver"],
    threads:
        ANN_DIA.get("threads", 8)
    resources:
        mem_mb = CFG_RES.get("diamond_mem_mb", 8192),
    conda:
        "../envs/annotation.yaml"
    log:
        "logs/02_annotation/diamond__{group}.log"
    shell:
        r"""
        bash workflow/scripts/15_diamond_annotate.sh \
            "{input.pos}" "{input.neg}" \
            "{input.db_uniprot}" "{input.db_ncbi}" "{input.db_conoserver}" \
            "{params.outdir}" "{params.evalue}" "{params.qcov}" \
            "{params.hsps}" "{params.mts_up}" "{params.mts_nc}" \
            "{params.mts_cs}" "{threads}" \
            > {log} 2>&1
        """


rule parse_besthits:
    input:
        outs = expand(
            f"{D_ANN}/{{{{group}}}}/annotate/blast/{{part}}__{{db}}.outfmt6",
            part=["positivo", "negativo"],
            db=["uniprot", "ncbi", "conoserver"]),
    output:
        besthits = f"{D_ANN}/{{group}}/annotate/besthits.tsv",
        audit    = f"{D_ANN}/{{group}}/annotate/blast_audit.tsv",
        flagged  = f"{D_ANN}/{{group}}/annotate/blast_flagged.tsv",
    params:
        blast_dir = lambda w: f"{D_ANN}/{w.group}/annotate/blast",
        tsv_dir   = lambda w: f"{D_ANN}/{w.group}/annotate/blast_tsv",
        qcov      = ANN_BH["min_qcov"],
        cov_delta = ANN_BH["cov_delta"],
        ratio     = ANN_BH["bitscore_ratio"],
    resources:
        mem_mb = CFG_RES.get("py_mem_mb", 4096),
    conda:
        "../envs/python_base.yaml"
    log:
        "logs/02_annotation/besthits__{group}.log"
    shell:
        r"""
        python workflow/scripts/16_parse_besthits.py \
            --blast-dir "{params.blast_dir}" \
            --out-besthits "{output.besthits}" \
            --out-tsv-dir "{params.tsv_dir}" \
            --out-audit "{output.audit}" --out-flagged "{output.flagged}" \
            --min-qcov "{params.qcov}" --cov-delta "{params.cov_delta}" \
            --bitscore-ratio "{params.ratio}" > {log} 2>&1
        """


# ============================================================================
# R7 — label resolution
# ============================================================================
rule annotate_labels:
    input:
        besthits = f"{D_ANN}/{{group}}/annotate/besthits.tsv",
        pos_xlsx = f"{D_ANN}/{{group}}/signalp/positivo_superfamilia.xlsx",
        neg_xlsx = f"{D_ANN}/{{group}}/signalp/negativo_superfamilia.xlsx",
        cs_fasta = _db("conoserver", "fasta"),
        ncbi_ft  = f"{DB_ROOT}/{ANN_DB['ncbi']['feature_table']}",
        up_xlsx  = f"{DB_ROOT}/{ANN_DB['uniprot']['labels_xlsx']}",
    output:
        pos    = f"{D_ANN}/{{group}}/annotate/positivo_annotated.xlsx",
        neg    = f"{D_ANN}/{{group}}/annotate/negativo_annotated.xlsx",
        labels = f"{D_ANN}/{{group}}/annotate/labels.tsv",
        audit  = f"{D_ANN}/{{group}}/annotate/annotate_audit.tsv",
    params:
        cs_field = ANN_DB["conoserver"]["header_field"],
        up_sheet = ANN_DB["uniprot"]["labels_sheet"],
        up_entry = ANN_DB["uniprot"]["labels_entry_col"],
        up_name  = ANN_DB["uniprot"]["labels_name_col"],
    resources:
        mem_mb = CFG_RES.get("py_mem_mb", 4096),
    conda:
        "../envs/python_base.yaml"
    log:
        "logs/02_annotation/labels__{group}.log"
    shell:
        r"""
        python workflow/scripts/17_annotate_labels.py \
            --besthits "{input.besthits}" \
            --pos-xlsx "{input.pos_xlsx}" --neg-xlsx "{input.neg_xlsx}" \
            --conoserver-fasta "{input.cs_fasta}" \
            --conoserver-field "{params.cs_field}" \
            --ncbi-feature-table "{input.ncbi_ft}" \
            --uniprot-xlsx "{input.up_xlsx}" \
            --uniprot-sheet "{params.up_sheet}" \
            --uniprot-entry-col "{params.up_entry}" \
            --uniprot-name-col "{params.up_name}" \
            --key-col "seq_id" \
            --out-pos-xlsx "{output.pos}" --out-neg-xlsx "{output.neg}" \
            --out-labels "{output.labels}" --out-audit "{output.audit}" \
            > {log} 2>&1
        """


# ============================================================================
# R8/R9 — normalisation and confidence tiers
# ============================================================================
rule normalize_final:
    input:
        pos = f"{D_ANN}/{{group}}/annotate/positivo_annotated.xlsx",
        neg = f"{D_ANN}/{{group}}/annotate/negativo_annotated.xlsx",
    output:
        pos_x      = f"{D_ANN}/{{group}}/classify/positivo_normalized.xlsx",
        neg_x      = f"{D_ANN}/{{group}}/classify/negativo_normalized.xlsx",
        normalized = f"{D_ANN}/{{group}}/classify/normalized.tsv",
        unresolved = f"{D_ANN}/{{group}}/classify/unresolved_labels.tsv",
        audit      = f"{D_ANN}/{{group}}/classify/normalize_audit.tsv",
    resources:
        mem_mb = CFG_RES.get("py_mem_mb", 4096),
    conda:
        "../envs/python_base.yaml"
    log:
        "logs/02_annotation/normalize__{group}.log"
    shell:
        r"""
        python workflow/scripts/18_normalize_final.py \
            --pos-xlsx "{input.pos}" --neg-xlsx "{input.neg}" \
            --out-pos-xlsx "{output.pos_x}" --out-neg-xlsx "{output.neg_x}" \
            --out-normalized "{output.normalized}" \
            --out-unresolved "{output.unresolved}" \
            --out-audit "{output.audit}" > {log} 2>&1
        """


rule classify_confidence:
    input:
        normalized = f"{D_ANN}/{{group}}/classify/normalized.tsv",
    output:
        xlsx    = f"{D_ANN}/{{group}}/classify/classification.xlsx",
        tsv     = f"{D_ANN}/{{group}}/classify/classification.tsv",
        summary = f"{D_ANN}/{{group}}/classify/tier_summary.tsv",
        audit   = f"{D_ANN}/{{group}}/classify/classify_audit.tsv",
    resources:
        mem_mb = CFG_RES.get("py_mem_mb", 4096),
    conda:
        "../envs/python_base.yaml"
    log:
        "logs/02_annotation/classify__{group}.log"
    shell:
        r"""
        python workflow/scripts/19_classify_confidence.py \
            --normalized "{input.normalized}" \
            --out-xlsx "{output.xlsx}" --out-tsv "{output.tsv}" \
            --out-summary "{output.summary}" --out-audit "{output.audit}" \
            > {log} 2>&1
        """


# ============================================================================
# R10 — manifest, invariants, gate
# ============================================================================
def _all_tier_summaries(wildcards=None):
    return [f"{D_ANN}/{g}/classify/tier_summary.tsv" for g in group_ids()]


rule annotation_manifest:
    input:
        groups    = GROUPS_VALIDATED,
        summaries = _all_tier_summaries,
    output:
        manifest = f"{D_ANN}/annotation_manifest.tsv",
        sha256   = f"{D_ANN}/annotation_manifest.tsv.sha256",
    params:
        root       = D_ANN,
        assemblers = ",".join(enabled_assemblers()),
    resources:
        mem_mb = CFG_RES.get("py_mem_mb", 4096),
    conda:
        "../envs/python_base.yaml"
    log:
        "logs/02_annotation/manifest.log"
    shell:
        r"""
        python workflow/scripts/20_aggregate_manifest.py \
            --groups "{input.groups}" --results-root "{params.root}" \
            --assemblers "{params.assemblers}" \
            --out-manifest "{output.manifest}" > {log} 2>&1
        """


rule annotation_invariants:
    input:
        manifest = f"{D_ANN}/annotation_manifest.tsv",
    output:
        report  = f"{D_QC}/annotation_invariants.txt",
        verdict = f"{D_QC}/annotation_verdict.txt",
    params:
        min_rows = ANN_INV.get("require_manifest_min_rows", 1),
        req_nr   = str(ANN_INV.get("require_nonempty_nr_pep", True)).lower(),
        req_sp   = str(ANN_INV.get("require_nonempty_signalp", False)).lower(),
        max_pct  = ANN_INV.get("max_conflicting_labels_pct", 25.0),
    conda:
        "../envs/python_base.yaml"
    log:
        "logs/02_annotation/invariants.log"
    shell:
        r"""
        python workflow/scripts/21_check_invariants.py \
            --manifest "{input.manifest}" \
            --out-report "{output.report}" --out-verdict "{output.verdict}" \
            --min-rows "{params.min_rows}" \
            --require-nonempty-nr-pep "{params.req_nr}" \
            --require-nonempty-signalp "{params.req_sp}" \
            --max-conflicting-pct "{params.max_pct}" > {log} 2>&1
        """


rule pipeline_gate:
    input:
        verdict = f"{D_QC}/annotation_verdict.txt",
        report  = f"{D_QC}/annotation_invariants.txt",
    output:
        ok = f"{D_QC}/pipeline_gate.OK",
    shell:
        r"""
        v=$(head -1 {input.verdict})
        if [ "$v" != "PASS" ]; then
            echo "Pipeline gate: verdict=$v — halting." >&2
            echo "See {input.report}" >&2
            exit 1
        fi
        touch {output.ok}
        """


# ---- Stage targets ----------------------------------------------------------
def stage02_targets():
    return [
        f"{D_ANN}/annotation_manifest.tsv",
        f"{D_QC}/annotation_invariants.txt",
        f"{D_QC}/pipeline_gate.OK",
    ]


localrules: annotation_preflight, annotation_preflight_gate, \
            annotation_manifest, annotation_invariants, pipeline_gate
