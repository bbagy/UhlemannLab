###############################################
# Go_MAGs_QC_V1.smk
#
# QC + host read depletion for MAG assembly input.
# Output names intentionally match the legacy GoQC.smk contract:
#   filtered_fastq/<sample>_R1_filtered.fastq.gz
#   filtered_fastq/<sample>_R2_filtered.fastq.gz
#   host_filtered_fastq/<sample>_R1_nohuman.fastq.gz
#   host_filtered_fastq/<sample>_R2_nohuman.fastq.gz
#   QC_summary.csv
###############################################

import glob
import gzip
import json
import os
import re
from pathlib import Path

# ---------- Config ----------
FASTQ_DIR = config.get("fastq_dir", config.get("Fastq_DIRS", "input_fastqs"))
OUTPUT_DIR = config.get("output_dir", config.get("output", "output"))
HOST_DB = config.get("host_db", config.get("hostDB", None))
THREADS = int(config.get("threads", 8))
PAIRED = int(config.get("paired", 2))

if HOST_DB is None:
    raise ValueError("[Go_MAGs_QC] Missing host DB. Use --config host_db=/path/to/bowtie2/index_prefix")
if PAIRED != 2:
    raise ValueError("[Go_MAGs_QC] MAGs QC currently expects paired-end reads: --config paired=2")

FILTERED_DIR = f"{OUTPUT_DIR}/filtered_fastq"
HOST_FILTERED_DIR = f"{OUTPUT_DIR}/host_filtered_fastq"
INTERMEDIATE_DIR = f"{OUTPUT_DIR}/intermediate"
LOG_DIR = f"{OUTPUT_DIR}/logs"

FASTQ_PATTERNS = [
    "*_L00[1-9]_R1_001.fastq.gz",
    "*_R1_001.fastq.gz",
    "*_R1.fastq.gz",
    "*.R1.fastq.gz",
]


def _sample_from_r1(path):
    base = os.path.basename(path)
    for pattern in [
        r"(.+)_L00[1-9]_R1_001\.fastq\.gz$",
        r"(.+)_R1_001\.fastq\.gz$",
        r"(.+)_R1\.fastq\.gz$",
        r"(.+)\.R1\.fastq\.gz$",
    ]:
        match = re.match(pattern, base)
        if match:
            return match.group(1)
    return None


def _discover_samples():
    r1s = []
    for pat in FASTQ_PATTERNS:
        r1s.extend(glob.glob(os.path.join(FASTQ_DIR, pat)))
    samples = sorted({sample for path in r1s if (sample := _sample_from_r1(path))})
    if not samples:
        raise ValueError(f"[Go_MAGs_QC] No paired-end R1 FASTQs found in {FASTQ_DIR}")
    return samples


SAMPLES = _discover_samples()


def _r1_files(wildcards):
    sample = wildcards.sample
    patterns = [
        f"{FASTQ_DIR}/{sample}_L00[1-9]_R1_001.fastq.gz",
        f"{FASTQ_DIR}/{sample}_R1_001.fastq.gz",
        f"{FASTQ_DIR}/{sample}_R1.fastq.gz",
        f"{FASTQ_DIR}/{sample}.R1.fastq.gz",
    ]
    files = []
    for pat in patterns:
        files.extend(glob.glob(pat))
    files = sorted(set(files))
    if not files:
        raise ValueError(f"[Go_MAGs_QC] R1 FASTQ missing for {sample}")
    return files


def _r2_files(wildcards):
    sample = wildcards.sample
    r1s = _r1_files(wildcards)
    r2s = []
    for r1 in r1s:
        candidates = [
            r1.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz"),
            r1.replace("_R1.fastq.gz", "_R2.fastq.gz"),
            r1.replace(".R1.fastq.gz", ".R2.fastq.gz"),
        ]
        found = [p for p in candidates if p != r1 and os.path.exists(p)]
        if not found:
            raise ValueError(f"[Go_MAGs_QC] R2 pair missing for {r1}")
        r2s.append(found[0])
    return r2s


def _count_reads_gz(path):
    lines = 0
    with gzip.open(path, "rt") as fh:
        for lines, _ in enumerate(fh, start=1):
            pass
    return lines // 4


rule all:
    input:
        f"{OUTPUT_DIR}/QC_summary.csv",
        expand(f"{HOST_FILTERED_DIR}/{{sample}}_R1_nohuman.fastq.gz", sample=SAMPLES),
        expand(f"{HOST_FILTERED_DIR}/{{sample}}_R2_nohuman.fastq.gz", sample=SAMPLES)


rule make_dirs:
    output:
        touch(f"{OUTPUT_DIR}/.dirs.ok")
    run:
        for d in [OUTPUT_DIR, FILTERED_DIR, HOST_FILTERED_DIR, INTERMEDIATE_DIR, LOG_DIR]:
            os.makedirs(d, exist_ok=True)
        Path(output[0]).touch()


rule fastp_qc:
    input:
        dirs=f"{OUTPUT_DIR}/.dirs.ok",
        r1=_r1_files,
        r2=_r2_files
    output:
        r1=f"{FILTERED_DIR}/{{sample}}_R1_filtered.fastq.gz",
        r2=f"{FILTERED_DIR}/{{sample}}_R2_filtered.fastq.gz",
        json=f"{FILTERED_DIR}/{{sample}}_fastp.json",
        html=f"{FILTERED_DIR}/{{sample}}_fastp.html"
    log:
        f"{LOG_DIR}/{{sample}}_fastp.log"
    params:
        tmp_r1=lambda wc: f"{INTERMEDIATE_DIR}/{wc.sample}.merged_R1.fastq.gz",
        tmp_r2=lambda wc: f"{INTERMEDIATE_DIR}/{wc.sample}.merged_R2.fastq.gz",
        n_lanes=lambda wc, input: len(input.r1),
        first_r1=lambda wc, input: input.r1[0],
        first_r2=lambda wc, input: input.r2[0]
    threads: THREADS
    shell:
        r"""
        set -euo pipefail

        if [ {params.n_lanes} -eq 1 ]; then
            r1_in={params.first_r1:q}
            r2_in={params.first_r2:q}
        else
            cat {input.r1:q} > {params.tmp_r1:q}
            cat {input.r2:q} > {params.tmp_r2:q}
            r1_in={params.tmp_r1:q}
            r2_in={params.tmp_r2:q}
        fi

        fastp \
            -i "$r1_in" \
            -I "$r2_in" \
            -o {output.r1:q} \
            -O {output.r2:q} \
            --detect_adapter_for_pe \
            --thread {threads} \
            --length_required 36 \
            --qualified_quality_phred 20 \
            --cut_front \
            --cut_tail \
            --html {output.html:q} \
            --json {output.json:q} \
            > {log:q} 2>&1
        """


rule bowtie2_host_filter:
    input:
        r1=f"{FILTERED_DIR}/{{sample}}_R1_filtered.fastq.gz",
        r2=f"{FILTERED_DIR}/{{sample}}_R2_filtered.fastq.gz"
    output:
        r1=f"{HOST_FILTERED_DIR}/{{sample}}_R1_nohuman.fastq.gz",
        r2=f"{HOST_FILTERED_DIR}/{{sample}}_R2_nohuman.fastq.gz"
    log:
        f"{LOG_DIR}/{{sample}}_host_filter.log"
    params:
        host_db=HOST_DB,
        tmp_r1=lambda wc: f"{INTERMEDIATE_DIR}/{wc.sample}.R1_nohuman.fastq",
        tmp_r2=lambda wc: f"{INTERMEDIATE_DIR}/{wc.sample}.R2_nohuman.fastq"
    threads: THREADS
    shell:
        r"""
        set -euo pipefail

        bowtie2 \
            -x {params.host_db:q} \
            -1 {input.r1:q} \
            -2 {input.r2:q} \
            --threads {threads} \
            --very-sensitive-local \
            --score-min L,2,0 \
            --no-mixed \
            --no-discordant \
            2> {log:q} \
        | samtools view -@ {threads} -b -f 12 -F 256 - \
            2>> {log:q} \
        | samtools sort -@ {threads} -n - \
            2>> {log:q} \
        | samtools fastq -@ {threads} \
            -1 {params.tmp_r1:q} \
            -2 {params.tmp_r2:q} \
            -0 /dev/null \
            -s /dev/null \
            -n - \
            2>> {log:q}

        pigz -p {threads} -c {params.tmp_r1:q} > {output.r1:q}
        pigz -p {threads} -c {params.tmp_r2:q} > {output.r2:q}
        rm -f {params.tmp_r1:q} {params.tmp_r2:q}
        """


rule qc_summary:
    input:
        jsons=expand(f"{FILTERED_DIR}/{{sample}}_fastp.json", sample=SAMPLES),
        host_r1=expand(f"{HOST_FILTERED_DIR}/{{sample}}_R1_nohuman.fastq.gz", sample=SAMPLES)
    output:
        summary=f"{OUTPUT_DIR}/QC_summary.csv"
    run:
        rows = ["Sample,Before Reads,After QC Reads,After Host Filter Reads"]
        for sample in SAMPLES:
            json_file = f"{FILTERED_DIR}/{sample}_fastp.json"
            host_r1 = f"{HOST_FILTERED_DIR}/{sample}_R1_nohuman.fastq.gz"
            with open(json_file) as fh:
                data = json.load(fh)
            before = data["summary"]["before_filtering"]["total_reads"]
            after_qc = data["summary"]["after_filtering"]["total_reads"]
            after_host_filter = _count_reads_gz(host_r1) if os.path.exists(host_r1) else 0
            rows.append(f"{sample},{before},{after_qc},{after_host_filter}")
        with open(output.summary, "w") as out:
            out.write("\n".join(rows) + "\n")
