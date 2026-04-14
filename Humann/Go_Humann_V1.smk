###############################################
# Go_Humann_V1.smk
#
# Snakemake rewrite of the HUMAnN3 execution block
# used in Uhlemann Lab shotgun workflows.
# Output contract is intentionally kept close to the
# shell history:
#   1_humann3_out/
#   2_humann3_final_out/
#   3_MetaPhlAn_bug_list/
#   4_kegg-orthology/
#   5_pathabundance_stratified_out/
#   humann3_log.txt
###############################################

import glob
import os
import re
import shutil
from pathlib import Path

FASTQ_DIR = config.get("fastq_dir", "host_filtered_fastq")
OUTPUT_DIR = config.get("output_dir", config.get("output", "output"))
NUCLEOTIDE_DB = config.get("nucleotide_db", "/media/uhlemann/core4/DB/humann_db/humann3/chocophlan")
PROTEIN_DB = config.get("protein_db", "/media/uhlemann/core4/DB/humann_db/humann3/uniref")
METAPHLAN_DB = config.get("metaphlan_db", "")
METAPHLAN_INDEX = config.get("metaphlan_index", "")
THREADS_PER_SAMPLE = int(config.get("humann_threads", config.get("threads", 4)))
KEEP_LOGS = int(config.get("keep_logs", 1))
RUN_GENE_NORM = str(config.get("run_gene_norm", "true")).lower() in ["1", "true", "yes", "y"]
RUN_PATH_SPLIT = str(config.get("run_path_split", "true")).lower() in ["1", "true", "yes", "y"]
RUN_PATHCOVERAGE_MERGE = str(config.get("run_pathcoverage_merge", "true")).lower() in ["1", "true", "yes", "y"]
RUN_MUSICC = str(config.get("run_musicc", "false")).lower() in ["1", "true", "yes", "y"]

RUN_DIR = f"{OUTPUT_DIR}/1_humann3_out"
FINAL_DIR = f"{OUTPUT_DIR}/2_humann3_final_out"
BUG_DIR = f"{OUTPUT_DIR}/3_MetaPhlAn_bug_list"
KO_DIR = f"{OUTPUT_DIR}/4_kegg-orthology"
KO_STRAT_DIR = f"{KO_DIR}/stratified_out"
PATH_STRAT_DIR = f"{OUTPUT_DIR}/5_pathabundance_stratified_out"
LOG_DIR = f"{OUTPUT_DIR}/6_logs"
INTERMEDIATE_DIR = f"{OUTPUT_DIR}/7_intermediate"

SCRIPT_DIR = os.path.join(workflow.basedir, "scripts")

FASTQ_PATTERNS = ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]


def _validate_runtime():
    missing_paths = []
    for label, path in [
        ("nucleotide_db", NUCLEOTIDE_DB),
        ("protein_db", PROTEIN_DB),
    ]:
        if not path or not os.path.exists(path):
            missing_paths.append(f"{label}={path}")

    if METAPHLAN_DB and not os.path.exists(METAPHLAN_DB):
        missing_paths.append(f"metaphlan_db={METAPHLAN_DB}")

    missing_bins = [exe for exe in ["humann", "metaphlan"] if shutil.which(exe) is None]

    if not missing_paths and not missing_bins:
        return

    msg = ["[Humann] Runtime validation failed."]
    if missing_paths:
        msg.append("Missing paths: " + ", ".join(missing_paths))
    if missing_bins:
        msg.append("Missing executables in PATH: " + ", ".join(missing_bins))
    msg.append(
        "If you are using the Docker workflow, launch through Go_Humann.sh so host paths are mounted "
        "to /fastq and /db inside the container. If you are running Snakemake directly, pass real host "
        "paths with --config nucleotide_db=/path/to/chocophlan protein_db=/path/to/uniref "
        "[metaphlan_db=/path/to/metaphlan_db]."
    )
    raise ValueError(" ".join(msg))


def _sample_from_fastq(path):
    base = os.path.basename(path)
    stem = re.sub(r"\.(fastq|fq)(\.gz)?$", "", base, flags=re.IGNORECASE)
    paired_patterns = [
        r"(.+)_R[12]_nohuman$",
        r"(.+)_R[12]$",
        r"(.+)\.R[12]$",
        r"(.+)_L00[1-9]_R[12]_001$",
        r"(.+)_R[12]_001$",
    ]
    for pattern in paired_patterns:
        match = re.match(pattern, stem, flags=re.IGNORECASE)
        if match:
            return match.group(1)
    return stem


def _discover_fastqs():
    fastqs = []
    for pat in FASTQ_PATTERNS:
        fastqs.extend(glob.glob(os.path.join(FASTQ_DIR, pat)))
    fastqs = sorted(set(fastqs))
    if not fastqs:
        raise ValueError(f"[Humann] No FASTQ files found in {FASTQ_DIR}")

    sample_map = {}
    for f in fastqs:
        sample = _sample_from_fastq(f)
        sample_map.setdefault(sample, []).append(f)

    for sample, files in sample_map.items():
        if len(files) == 1:
            continue
        if len(files) != 2:
            raise ValueError(f"[Humann] Expected 1 or 2 FASTQs for {sample}, found {len(files)}: {files}")

        names = [os.path.basename(f) for f in sorted(files)]
        is_pair = (
            any(re.search(r"_R1(_nohuman)?(\.|$)|\.R1(\.|$)|_R1_001", name, flags=re.IGNORECASE) for name in names)
            and any(re.search(r"_R2(_nohuman)?(\.|$)|\.R2(\.|$)|_R2_001", name, flags=re.IGNORECASE) for name in names)
        )
        if not is_pair:
            raise ValueError(f"[Humann] Multiple FASTQs for {sample} do not look like an R1/R2 pair: {files}")

        sample_map[sample] = sorted(
            files,
            key=lambda x: 1 if re.search(r"_R1(_nohuman)?(\.|$)|\.R1(\.|$)|_R1_001", os.path.basename(x), flags=re.IGNORECASE) else 2
        )
    return sample_map


SAMPLE_FASTQS = _discover_fastqs()
SAMPLES = sorted(SAMPLE_FASTQS.keys())
_validate_runtime()


def _fastq(wildcards):
    return SAMPLE_FASTQS[wildcards.sample]


def _humann_input(wildcards):
    return f"{INTERMEDIATE_DIR}/{wildcards.sample}.fastq"


ALL_TARGETS = [
    f"{FINAL_DIR}/merged_genefamilies.txt",
    f"{FINAL_DIR}/merged_pathabundance.txt",
        f"{OUTPUT_DIR}/humann3_log.txt",
    *expand(f"{RUN_DIR}/{{sample}}_genefamilies.tsv", sample=SAMPLES),
    *expand(f"{RUN_DIR}/{{sample}}_pathabundance.tsv", sample=SAMPLES),
    *expand(f"{BUG_DIR}/{{sample}}_metaphlan_bugs_list.tsv", sample=SAMPLES),
    *expand(f"{INTERMEDIATE_DIR}/{{sample}}.input.done", sample=SAMPLES),
]

if RUN_PATHCOVERAGE_MERGE:
    ALL_TARGETS.extend([
        f"{FINAL_DIR}/merged_pathcoverage.txt",
        *expand(f"{RUN_DIR}/{{sample}}_pathcoverage.tsv", sample=SAMPLES),
    ])

if RUN_GENE_NORM:
    ALL_TARGETS.extend([
        f"{KO_DIR}/merged_genefamilies_cpm.txt",
        f"{KO_DIR}/merged_genefamilies_uniref90_rxn_cpm.txt",
        f"{KO_DIR}/merged_genefamilies_uniref90_rxn_musicc.txt" if RUN_MUSICC else f"{KO_DIR}/.musicc.skip",
        f"{KO_DIR}/merged_genefamilies_uniref90_rxn_kegg-orthology_cpm.txt",
        f"{KO_STRAT_DIR}/merged_genefamilies_uniref90_rxn_kegg-orthology_cpm_unstratified.txt",
        f"{KO_STRAT_DIR}/merged_genefamilies_uniref90_rxn_kegg-orthology_cpm_unstratified_filtered.txt",
    ])

if RUN_PATH_SPLIT:
    ALL_TARGETS.extend([
        f"{PATH_STRAT_DIR}/merged_pathabundance_unstratified.txt",
        f"{PATH_STRAT_DIR}/.split.done",
    ])


rule all:
    input:
        ALL_TARGETS


rule make_dirs:
    output:
        touch(f"{OUTPUT_DIR}/.dirs.ok")
    run:
        for d in [
            OUTPUT_DIR, RUN_DIR, FINAL_DIR, BUG_DIR, KO_DIR, KO_STRAT_DIR, PATH_STRAT_DIR, LOG_DIR, INTERMEDIATE_DIR
        ]:
            os.makedirs(d, exist_ok=True)
        Path(output[0]).touch()


rule prepare_input:
    input:
        dirs=f"{OUTPUT_DIR}/.dirs.ok",
        fastq=_fastq
    output:
        merged=temp(f"{INTERMEDIATE_DIR}/{{sample}}.fastq"),
        done=touch(f"{INTERMEDIATE_DIR}/{{sample}}.input.done")
    params:
        r1=lambda wc, input: input.fastq[0] if isinstance(input.fastq, list) else input.fastq,
        r2=lambda wc, input: input.fastq[1] if isinstance(input.fastq, list) and len(input.fastq) > 1 else ""
    shell:
        r"""
        set -euo pipefail

        if [ -n "{params.r2}" ]; then
            zcat {input.fastq:q} > {output.merged:q}
        else
            zcat {params.r1:q} > {output.merged:q}
        fi
        """


rule run_humann:
    input:
        dirs=f"{OUTPUT_DIR}/.dirs.ok",
        prepared=f"{INTERMEDIATE_DIR}/{{sample}}.input.done",
        fastq=_humann_input
    output:
        genefamilies=f"{RUN_DIR}/{{sample}}_genefamilies.tsv",
        pathabundance=f"{RUN_DIR}/{{sample}}_pathabundance.tsv",
        pathcoverage=f"{RUN_DIR}/{{sample}}_pathcoverage.tsv",
        bug_list=f"{BUG_DIR}/{{sample}}_metaphlan_bugs_list.tsv",
        done=touch(f"{RUN_DIR}/{{sample}}.humann.done")
    log:
        f"{LOG_DIR}/{{sample}}.humann.log"
    params:
        outdir=RUN_DIR,
        nucleotide_db=NUCLEOTIDE_DB,
        protein_db=PROTEIN_DB,
        metaphlan_db=METAPHLAN_DB,
        metaphlan_index=METAPHLAN_INDEX
    threads: THREADS_PER_SAMPLE
    shell:
        r"""
        set -euo pipefail

        echo "[run_humann] start $(date)" > {log:q}
        echo "[run_humann] sample={wildcards.sample}" >> {log:q}
        echo "[run_humann] input={input.fastq}" >> {log:q}
        echo "[run_humann] output_dir={params.outdir}" >> {log:q}
        echo "[run_humann] nucleotide_db={params.nucleotide_db}" >> {log:q}
        echo "[run_humann] protein_db={params.protein_db}" >> {log:q}
        echo "[run_humann] metaphlan_db={params.metaphlan_db}" >> {log:q}
        echo "[run_humann] metaphlan_index={params.metaphlan_index}" >> {log:q}
        echo "[run_humann] humann=$(command -v humann || true)" >> {log:q}
        echo "[run_humann] metaphlan=$(command -v metaphlan || true)" >> {log:q}
        humann --version >> {log:q} 2>&1 || true
        metaphlan --version >> {log:q} 2>&1 || true

        humann_args=(
            --threads {threads}
            --input {input.fastq:q}
            --output {params.outdir:q}
            --nucleotide-database {params.nucleotide_db:q}
            --protein-database {params.protein_db:q}
        )

        metaphlan_args=()
        if [ -n "{params.metaphlan_db}" ]; then
            metaphlan_args+=(--bowtie2db "{params.metaphlan_db}")
        fi
        if [ -n "{params.metaphlan_index}" ]; then
            metaphlan_args+=(--index "{params.metaphlan_index}")
        fi
        if [ "${{#metaphlan_args[@]}}" -gt 0 ]; then
            humann_args+=(--metaphlan-options "${{metaphlan_args[*]}}")
        fi

        printf '[run_humann] cmd: humann' >> {log:q}
        printf ' %q' "${{humann_args[@]}}" >> {log:q}
        printf '\n' >> {log:q}

        humann "${{humann_args[@]}}" >> {log:q} 2>&1

        temp_dir="{params.outdir}/{wildcards.sample}_humann_temp"
        bugs_src="$temp_dir/{wildcards.sample}_metaphlan_bugs_list.tsv"
        if [ -f "$bugs_src" ]; then
            mv "$bugs_src" {output.bug_list:q}
        else
            : > {output.bug_list:q}
        fi
        rm -rf "$temp_dir"
        """


rule merge_genefamilies:
    input:
        expand(f"{RUN_DIR}/{{sample}}_genefamilies.tsv", sample=SAMPLES)
    output:
        f"{FINAL_DIR}/merged_genefamilies.txt"
    log:
        f"{LOG_DIR}/merge_genefamilies.log"
    shell:
        r"""
        set -euo pipefail
        humann_join_tables \
            -i {RUN_DIR:q} \
            -o {output:q} \
            --file_name genefamilies \
            > {log:q} 2>&1
        """


rule merge_pathabundance:
    input:
        expand(f"{RUN_DIR}/{{sample}}_pathabundance.tsv", sample=SAMPLES)
    output:
        f"{FINAL_DIR}/merged_pathabundance.txt"
    log:
        f"{LOG_DIR}/merge_pathabundance.log"
    shell:
        r"""
        set -euo pipefail
        humann_join_tables \
            -i {RUN_DIR:q} \
            -o {output:q} \
            --file_name pathabundance \
            > {log:q} 2>&1
        """


rule merge_pathcoverage:
    input:
        expand(f"{RUN_DIR}/{{sample}}_pathcoverage.tsv", sample=SAMPLES)
    output:
        f"{FINAL_DIR}/merged_pathcoverage.txt"
    log:
        f"{LOG_DIR}/merge_pathcoverage.log"
    shell:
        r"""
        set -euo pipefail
        humann_join_tables \
            -i {RUN_DIR:q} \
            -o {output:q} \
            --file_name pathcoverage \
            > {log:q} 2>&1
        """


rule gene_family_normalization:
    input:
        merged=f"{FINAL_DIR}/merged_genefamilies.txt"
    output:
        cpm=f"{KO_DIR}/merged_genefamilies_cpm.txt",
        regroup=f"{KO_DIR}/merged_genefamilies_uniref90_rxn_cpm.txt",
        musicc=f"{KO_DIR}/merged_genefamilies_uniref90_rxn_musicc.txt" if RUN_MUSICC else f"{KO_DIR}/.musicc.skip",
        renamed=f"{KO_DIR}/merged_genefamilies_uniref90_rxn_kegg-orthology_cpm.txt",
        unstrat=f"{KO_STRAT_DIR}/merged_genefamilies_uniref90_rxn_kegg-orthology_cpm_unstratified.txt",
        filtered=f"{KO_STRAT_DIR}/merged_genefamilies_uniref90_rxn_kegg-orthology_cpm_unstratified_filtered.txt"
    params:
        run_musicc=RUN_MUSICC
    log:
        f"{LOG_DIR}/gene_family_normalization.log"
    shell:
        r"""
        set -euo pipefail
        humann_renorm_table \
            --input {input.merged:q} \
            --output {output.cpm:q} \
            --units cpm \
            --update-snames \
            > {log:q} 2>&1

        humann_regroup_table \
            --input {output.cpm:q} \
            --output {output.regroup:q} \
            --groups uniref90_rxn \
            >> {log:q} 2>&1

        if [ "{params.run_musicc}" = "True" ]; then
            run_musicc.py {output.regroup:q} \
                -n -c use_generic -v \
                -o {output.musicc:q} \
                >> {log:q} 2>&1
            rename_input={output.musicc:q}
        else
            : > {output.musicc:q}
            rename_input={output.regroup:q}
        fi

        humann_rename_table \
            --input "$rename_input" \
            --output {output.renamed:q} \
            --names kegg-orthology \
            >> {log:q} 2>&1

        humann_split_stratified_table \
            --input {output.renamed:q} \
            --output {KO_STRAT_DIR:q} \
            >> {log:q} 2>&1

        head -n 1 {output.unstrat:q} > {output.filtered:q}
        grep ":" {output.unstrat:q} >> {output.filtered:q} || true
        """


rule split_pathabundance:
    input:
        merged=f"{FINAL_DIR}/merged_pathabundance.txt"
    output:
        unstrat=f"{PATH_STRAT_DIR}/merged_pathabundance_unstratified.txt",
        done=touch(f"{PATH_STRAT_DIR}/.split.done")
    log:
        f"{LOG_DIR}/split_pathabundance.log"
    shell:
        r"""
        set -euo pipefail
        humann_split_stratified_table \
            --input {input.merged:q} \
            --output {PATH_STRAT_DIR:q} \
            > {log:q} 2>&1
        """


rule master_log:
    input:
        expand(f"{LOG_DIR}/{{sample}}.humann.log", sample=SAMPLES)
    output:
        f"{OUTPUT_DIR}/humann3_log.txt"
    params:
        keep_logs=KEEP_LOGS
    shell:
        r"""
        set -euo pipefail
        python {SCRIPT_DIR}/humann_masterlog.py \
            --output {output:q} \
            --keep-logs {params.keep_logs} \
            {input:q}
        """
