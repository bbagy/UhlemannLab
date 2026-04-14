###############################################
# Go_KBracken_V1.smk
#
# Snakemake rewrite of Gobracken2_V4.pl.
# Output contract is intentionally kept close to the Perl wrapper:
#   1_out, 2_report, 3_classified, 4_unclassified,
#   5_mpa_report, 6_bracken_out, 7_bracken_mpa,
#   kraken2_mpa.txt, bracken_mpa.txt, bracken_mpa_filled.txt,
#   kraken2_log.txt
###############################################

import glob
import os
import re
from pathlib import Path

FASTQ_DIR = config.get("fastq_dir", "input_fastqs")
OUTPUT_DIR = config.get("output_dir", config.get("output", "output"))
DB = config.get("db", "/media/uhlemann/core4/DB/kraken2DB/k2_pluspfp_16gb_20241228")

DELIMITER = config.get("delimiter", "_")
if DELIMITER not in ["_", "."]:
    raise ValueError('config delimiter must be "_" or "."')

KRAKEN_THREADS = int(config.get("kraken_threads", config.get("threads", 4)))
BRACKEN_READ_LEN = int(config.get("bracken_read_len", 100))
BRACKEN_LEVEL = config.get("bracken_level", "S")
BRACKEN_THRESHOLD = int(config.get("bracken_threshold", 10))
RUN_BRACKEN = str(config.get("run_bracken", "true")).lower() in ["1", "true", "yes", "y"]
KEEP_LOGS = int(config.get("keep_logs", 1))

OUT_RAW_DIR = f"{OUTPUT_DIR}/1_out"
REPORT_DIR = f"{OUTPUT_DIR}/2_report"
CLASSIFIED_DIR = f"{OUTPUT_DIR}/3_classified"
UNCLASSIFIED_DIR = f"{OUTPUT_DIR}/4_unclassified"
MPA_REPORT_DIR = f"{OUTPUT_DIR}/5_mpa_report"
BRACKEN_DIR = f"{OUTPUT_DIR}/6_bracken_out"
BRACKEN_MPA_DIR = f"{OUTPUT_DIR}/7_bracken_mpa"
LOG_DIR = f"{OUTPUT_DIR}/logs"

SCRIPT_DIR = os.path.join(workflow.basedir, "scripts")

FASTQ_PATTERNS = ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]
R1_PAT = re.compile(r"(^|[_.])(R1|1)([_.]|$)|forward", re.IGNORECASE)
R2_PAT = re.compile(r"(^|[_.])(R2|2)([_.]|$)|reverse", re.IGNORECASE)
READ_TOKEN_PAT = re.compile(r"(^|[_.])(R[12])([_.]|$)|forward|reverse", re.IGNORECASE)


def _sample_from_fastq(path):
    base = os.path.basename(path)
    stem = re.sub(r"\.(fastq|fq)(\.gz)?$", "", base, flags=re.IGNORECASE)
    m = re.match(r"(.+?)[_.]R[12]([_.-].*)?$", stem, flags=re.IGNORECASE)
    if m:
        return m.group(1)
    if re.search("forward|reverse", stem, flags=re.IGNORECASE):
        return re.sub(r"([_.-])?(forward|reverse)([_.-].*)?$", "", stem, flags=re.IGNORECASE)
    return stem


def _discover_fastqs():
    fastqs = []
    for pat in FASTQ_PATTERNS:
        fastqs.extend(glob.glob(os.path.join(FASTQ_DIR, pat)))
    fastqs = sorted(set(fastqs))
    if not fastqs:
        raise ValueError(f"[KBracken] No FASTQ files found in {FASTQ_DIR}")

    sample_map = {}
    for f in fastqs:
        sample_map.setdefault(_sample_from_fastq(f), []).append(f)

    se = 0
    pe = 0
    for sample, files in list(sample_map.items()):
        if len(files) == 1:
            se += 1
            continue
        if not all(READ_TOKEN_PAT.search(os.path.basename(f)) for f in files):
            se += len(files)
            for f in files:
                sample_map[_sample_from_fastq(f)] = [f]
            del sample_map[sample]
            continue
        if len(files) == 2:
            pe += 1
            f0 = os.path.basename(files[0])
            f1 = os.path.basename(files[1])
            if R1_PAT.search(f0) and R2_PAT.search(f1):
                continue
            if R1_PAT.search(f1) and R2_PAT.search(f0):
                sample_map[sample] = [files[1], files[0]]
                continue
            raise ValueError(f"[KBracken] Paired files for {sample} do not look like R1/R2: {files}")
        raise ValueError(f"[KBracken] {len(files)} FASTQs found for sample {sample}: {files}")

    if se > 0 and pe > 0:
        raise ValueError(f"[KBracken] Mixed SE/PE input detected: SE={se}, PE={pe}. Run separately.")

    paired_config = config.get("paired", None)
    if paired_config is None:
        paired = pe > 0
    else:
        paired = str(paired_config).lower() in ["1", "true", "yes", "y", "paired", "2"]

    if paired and se > 0:
        raise ValueError(f"[KBracken] paired=true but {se} SE sample(s) detected.")
    if (not paired) and pe > 0:
        raise ValueError(f"[KBracken] paired=false but {pe} PE sample(s) detected.")

    return sample_map, paired


SAMPLE_FASTQS, PAIRED = _discover_fastqs()
SAMPLES = sorted(SAMPLE_FASTQS.keys())
PAIRED_FLAG = "--paired" if PAIRED else ""


def _fastqs(wildcards):
    return SAMPLE_FASTQS[wildcards.sample]


ALL_TARGETS = [
    f"{OUTPUT_DIR}/kraken2_mpa.txt",
    f"{OUTPUT_DIR}/kraken2_log.txt",
    *expand(f"{OUT_RAW_DIR}/{{sample}}_out.txt", sample=SAMPLES),
    *expand(f"{REPORT_DIR}/{{sample}}_report.txt", sample=SAMPLES),
    *expand(f"{MPA_REPORT_DIR}/{{sample}}_mpa.txt", sample=SAMPLES),
]

if RUN_BRACKEN:
    ALL_TARGETS.extend([
        f"{OUTPUT_DIR}/bracken_mpa.txt",
        f"{OUTPUT_DIR}/bracken_mpa_filled.txt",
        *expand(f"{BRACKEN_DIR}/{{sample}}_bracken.txt", sample=SAMPLES),
        *expand(f"{BRACKEN_MPA_DIR}/{{sample}}_bracken_mpa.txt", sample=SAMPLES),
    ])


rule all:
    input:
        ALL_TARGETS


rule make_dirs:
    output:
        touch(f"{OUTPUT_DIR}/.dirs.ok")
    run:
        for d in [
            OUTPUT_DIR, OUT_RAW_DIR, REPORT_DIR, CLASSIFIED_DIR, UNCLASSIFIED_DIR,
            MPA_REPORT_DIR, BRACKEN_DIR, BRACKEN_MPA_DIR, LOG_DIR
        ]:
            os.makedirs(d, exist_ok=True)
        Path(output[0]).touch()


rule kraken2_default:
    input:
        dirs=f"{OUTPUT_DIR}/.dirs.ok",
        fastqs=_fastqs
    output:
        out=f"{OUT_RAW_DIR}/{{sample}}_out.txt",
        report=f"{REPORT_DIR}/{{sample}}_report.txt",
        log=f"{OUT_RAW_DIR}/{{sample}}_out.log",
        done=touch(f"{OUT_RAW_DIR}/{{sample}}.kraken2.done")
    params:
        db=DB,
        paired_flag=PAIRED_FLAG,
        class_out=lambda wc: f"{CLASSIFIED_DIR}/{wc.sample}_classified#.fastq" if PAIRED else f"{CLASSIFIED_DIR}/{wc.sample}_classified.fastq",
        unclass_out=lambda wc: f"{UNCLASSIFIED_DIR}/{wc.sample}_unclassified#.fastq" if PAIRED else f"{UNCLASSIFIED_DIR}/{wc.sample}_unclassified.fastq",
        class_glob=lambda wc: f"{CLASSIFIED_DIR}/{wc.sample}_classified*.fastq",
        unclass_glob=lambda wc: f"{UNCLASSIFIED_DIR}/{wc.sample}_unclassified*.fastq"
    threads: KRAKEN_THREADS
    shell:
        r"""
        set -euo pipefail

        kraken2 {params.paired_flag} \
            --db {params.db:q} \
            --threads {threads} \
            --classified-out {params.class_out:q} \
            --unclassified-out {params.unclass_out:q} \
            --output {output.out:q} \
            --report {output.report:q} \
            {input.fastqs:q} \
            2> {output.log:q}

        for f in {params.class_glob} {params.unclass_glob}; do
            [ -e "$f" ] && pigz -f "$f"
        done
        """


rule kraken2_mpa:
    input:
        dirs=f"{OUTPUT_DIR}/.dirs.ok",
        fastqs=_fastqs
    output:
        mpa=f"{MPA_REPORT_DIR}/{{sample}}_mpa.txt"
    log:
        f"{LOG_DIR}/{{sample}}_kraken2_mpa.log"
    params:
        db=DB,
        paired_flag=PAIRED_FLAG
    threads: KRAKEN_THREADS
    shell:
        r"""
        set -euo pipefail
        kraken2 {params.paired_flag} \
            --db {params.db:q} \
            --threads {threads} \
            --use-mpa-style \
            --report {output.mpa:q} \
            --output /dev/null \
            {input.fastqs:q} \
            > /dev/null 2> {log:q}
        """


rule bracken:
    input:
        report=f"{REPORT_DIR}/{{sample}}_report.txt"
    output:
        bracken=f"{BRACKEN_DIR}/{{sample}}_bracken.txt"
    log:
        f"{LOG_DIR}/{{sample}}_bracken.log"
    params:
        db=DB,
        read_len=BRACKEN_READ_LEN,
        level=BRACKEN_LEVEL,
        threshold=BRACKEN_THRESHOLD
    shell:
        r"""
        set -euo pipefail
        bracken \
            -d {params.db:q} \
            -i {input.report:q} \
            -o {output.bracken:q} \
            -r {params.read_len} \
            -l {params.level:q} \
            -t {params.threshold} \
            > {log:q} 2>&1
        """


rule bracken_to_mpa:
    input:
        bracken=f"{BRACKEN_DIR}/{{sample}}_bracken.txt"
    output:
        mpa=f"{BRACKEN_MPA_DIR}/{{sample}}_bracken_mpa.txt"
    shell:
        r"""
        set -euo pipefail
        python {SCRIPT_DIR}/bracken_to_mpa.py \
            --input {input.bracken:q} \
            --output {output.mpa:q}
        """


rule merge_kraken_mpa:
    input:
        expand(f"{MPA_REPORT_DIR}/{{sample}}_mpa.txt", sample=SAMPLES)
    output:
        f"{OUTPUT_DIR}/kraken2_mpa.txt"
    shell:
        r"""
        set -euo pipefail
        python {SCRIPT_DIR}/merge_mpa_tables.py \
            --output {output:q} \
            {input:q}
        """


rule merge_bracken_mpa:
    input:
        expand(f"{BRACKEN_MPA_DIR}/{{sample}}_bracken_mpa.txt", sample=SAMPLES)
    output:
        f"{OUTPUT_DIR}/bracken_mpa.txt"
    shell:
        r"""
        set -euo pipefail
        python {SCRIPT_DIR}/merge_mpa_tables.py \
            --output {output:q} \
            {input:q}
        """


rule fill_bracken_taxonomy:
    input:
        kraken=f"{OUTPUT_DIR}/kraken2_mpa.txt",
        bracken=f"{OUTPUT_DIR}/bracken_mpa.txt"
    output:
        f"{OUTPUT_DIR}/bracken_mpa_filled.txt"
    shell:
        r"""
        set -euo pipefail
        python {SCRIPT_DIR}/fill_bracken_taxonomy.py \
            --kraken-mpa {input.kraken:q} \
            --bracken-mpa {input.bracken:q} \
            --output {output:q}
        """


rule master_log:
    input:
        expand(f"{OUT_RAW_DIR}/{{sample}}_out.log", sample=SAMPLES)
    output:
        f"{OUTPUT_DIR}/kraken2_log.txt"
    params:
        keep_logs=KEEP_LOGS
    shell:
        r"""
        set -euo pipefail
        python {SCRIPT_DIR}/kraken_masterlog.py \
            --output {output:q} \
            --keep-logs {params.keep_logs} \
            {input:q}
        """
