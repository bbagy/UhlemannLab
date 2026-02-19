###############################################
# Go_Pfsnake.smk  (PF WGS strain/variant typing)
#
# Modules:
# 1) fastp QC (flat output)
# 2) MultiQC
# 3) bwa mem -> sorted bam + index
# 4) Picard MarkDuplicates (optional, default ON) + BAM metrics
# 5) Variant calling (GATK HaplotypeCaller, haploid)
# 6) SNP summary (bcftools stats)
# 7) CNV calling:
#    - window-based depth CNV (log2 ratio)
#    - gene-based CNV (requires PlasmoDB GFF)
# 8) Mixed infection metrics:
#    - AF spectrum TSV
#    - heterozygosity proxy (haploid)
# 9) Recombination exploratory (subtelomeric):
#    - soft-clipped counts
#    - discordant pair counts (PE only)
# 10) Drug resistance gene summary:
#    - pfcrt / pfmdr1 / kelch13 variants extracted
# 11) Final reports (XLSX):
#    - YYYYMMDD_PF_QC_report.xlsx
#    - YYYYMMDD_PF_variant_summary.xlsx  (includes extra sheets)
###############################################

import os, re, glob, datetime, json
from pathlib import Path

# ----------------------------
# Required config keys
# ----------------------------
FASTQ_DIR  = config["Fastq_DIRS"]
OUTPUT_DIR = config["output"]
REF_FASTA  = config["hostDB"]
THREADS    = int(config.get("threads", 4))
PAIRED     = int(config.get("paired", 2))  # 2=PE, 1=SE

RUN_TAG = config.get("run_id", None) or datetime.datetime.now().strftime("%Y%m%d")

# ----------------------------
# Optional configs (CNV / GFF / Recomb / Drug)
# ----------------------------
# CNV window size (bp)
CNV_WINDOW = int(config.get("cnv_window", 1000))
# Subtelomeric region size (bp) for exploratory recomb
SUBTEL_SIZE = int(config.get("subtel_size", 50000))
# Use Picard MarkDuplicates (default ON)
USE_MARKDUP = int(config.get("use_markdup", 1))  # 1=ON, 0=OFF
# PlasmoDB GFF (for gene-level CNV + drug gene coordinate)
GFF_FILE = config.get("gff", None)  # e.g. /path/PlasmoDB-36_Pfalciparum3D7.gff

# Drug gene IDs (defaults are common 3D7 IDs; override if your GFF uses different IDs)
# You can also set these in config: drug_genes={"pfcrt":"PF3D7_0709000",...}
DRUG_GENES = config.get("drug_genes", {
    "pfcrt":   "PF3D7_0709000",
    "pfmdr1":  "PF3D7_0523000",
    "kelch13": "PF3D7_1343700",
})

# ----------------------------
# Output structure
# ----------------------------
QC_DIR      = f"{OUTPUT_DIR}/1_QC"
FASTP_DIR   = f"{QC_DIR}/fastp"
MULTIQC_DIR = f"{QC_DIR}/multiqc"

BWA_DIR     = f"{OUTPUT_DIR}/2_bwa"
BAMQC_DIR   = f"{OUTPUT_DIR}/3_bam_qc"

GATK_DIR    = f"{OUTPUT_DIR}/4_gatk_out"

SUM_DIR     = f"{OUTPUT_DIR}/5_variant_summary"
SUM_GATK    = f"{SUM_DIR}/1_gatk_summary"

CNV_DIR     = f"{OUTPUT_DIR}/6_cnv"
CNV_WIN_DIR = f"{CNV_DIR}/window"
CNV_GENE_DIR= f"{CNV_DIR}/gene"

MIX_DIR     = f"{OUTPUT_DIR}/7_mixed_infection"

RECOMB_DIR  = f"{OUTPUT_DIR}/8_recomb_explore"

DRUG_DIR    = f"{OUTPUT_DIR}/9_drug_resistance"

LOG_DIR     = f"{OUTPUT_DIR}/logs"
BENCH_DIR   = f"{OUTPUT_DIR}/benchmarks"
TMP_DIR     = f"{OUTPUT_DIR}/.tmp"

FINAL_QC_XLSX  = f"{OUTPUT_DIR}/{RUN_TAG}_PF_QC_report.xlsx"
FINAL_VAR_XLSX = f"{OUTPUT_DIR}/{RUN_TAG}_PF_variant_summary.xlsx"

# ----------------------------
# Sample discovery
# ----------------------------
def _find_samples(d):
    pats = ["*_R1.fastq.gz", "*_R1_001.fastq.gz", "*.R1.fastq.gz"]
    r1s = []
    for p in pats:
        r1s += sorted(glob.glob(os.path.join(d, p)))

    samples = set()
    for r1 in r1s:
        b = os.path.basename(r1)
        m = re.match(r"(.+?)(?:[_\.]R1(?:_001)?)\.fastq\.gz$", b)
        if m:
            samples.add(m.group(1))
    return sorted(samples)

SAMPLES = _find_samples(FASTQ_DIR)
if not SAMPLES:
    raise ValueError(f"[FATAL] No FASTQ samples found in {FASTQ_DIR}")

def _r1(sample):
    cands = [
        f"{FASTQ_DIR}/{sample}_R1_001.fastq.gz",
        f"{FASTQ_DIR}/{sample}_R1.fastq.gz",
        f"{FASTQ_DIR}/{sample}.R1.fastq.gz",
    ]
    for p in cands:
        if os.path.exists(p):
            return p
    raise ValueError(f"[FATAL] R1 missing for {sample}")

def _r2(sample):
    cands = [
        f"{FASTQ_DIR}/{sample}_R2_001.fastq.gz",
        f"{FASTQ_DIR}/{sample}_R2.fastq.gz",
        f"{FASTQ_DIR}/{sample}.R2.fastq.gz",
    ]
    for p in cands:
        if os.path.exists(p):
            return p
    raise ValueError(f"[FATAL] R2 missing for {sample}")

# ----------------------------
# Rule all (final targets)
# ----------------------------
rule all:
    input:
        # QC
        FINAL_QC_XLSX,

        # Variant summary
        FINAL_VAR_XLSX,

        # SNP tables
        f"{SUM_DIR}/all_samples.snps.long.tsv",
        f"{SUM_DIR}/all_samples.snps.wide.tsv",

        # Per-sample SNP tables (intermediate but keep for traceability)
        expand(
            f"{SUM_DIR}/snp_table/{{sample}}.snps.tsv",
            sample=SAMPLES
        ),

        # GATK stats
        expand(
            f"{SUM_GATK}/{{sample}}.bcftools.stats.txt",
            sample=SAMPLES
        ),

        # CNV (window)
        expand(
            f"{CNV_WIN_DIR}/{{sample}}.cnv.window.tsv",
            sample=SAMPLES
        ),

        # CNV (gene, optional)
        expand(
            f"{CNV_GENE_DIR}/{{sample}}.cnv.gene.tsv",
            sample=SAMPLES
        ) if GFF_FILE else [],

        # Mixed infection
        expand(
            f"{MIX_DIR}/{{sample}}.mixed.metrics.tsv",
            sample=SAMPLES
        ),

        # Recombination explore
        expand(
            f"{RECOMB_DIR}/{{sample}}.recomb_explore.tsv",
            sample=SAMPLES
        ),

        # Drug resistance (optional)
        expand(
            f"{DRUG_DIR}/{{sample}}.drug_variants.tsv",
            sample=SAMPLES
        ) if GFF_FILE else []

# ----------------------------
# Directory creation
# ----------------------------
rule make_dirs:
    output:
        touch(f"{OUTPUT_DIR}/.dirs.ok")
    run:
        for d in [
            OUTPUT_DIR, QC_DIR, FASTP_DIR, MULTIQC_DIR,
            BWA_DIR, BAMQC_DIR, GATK_DIR,
            SUM_DIR, SUM_GATK,
            CNV_DIR, CNV_WIN_DIR, CNV_GENE_DIR,
            MIX_DIR, RECOMB_DIR, DRUG_DIR,
            LOG_DIR, BENCH_DIR, TMP_DIR,
            f"{LOG_DIR}/fastp", f"{LOG_DIR}/bwa", f"{LOG_DIR}/picard", f"{LOG_DIR}/gatk",
            f"{LOG_DIR}/bcftools", f"{LOG_DIR}/cnv", f"{LOG_DIR}/mixed", f"{LOG_DIR}/recomb", f"{LOG_DIR}/drug"
        ]:
            os.makedirs(d, exist_ok=True)
        Path(output[0]).touch()

# ----------------------------
# Reference indexing
# ----------------------------
rule ref_index:
    input:
        ref=REF_FASTA
    output:
        bwa_ok = REF_FASTA + ".bwa.ok",
        fai    = REF_FASTA + ".fai",
        dict   = REF_FASTA.rsplit(".", 1)[0] + ".dict",
    log:
        f"{LOG_DIR}/ref_index.log"
    threads: 2
    shell:
        r"""
        set -euo pipefail

        # BWA index
        bwa index {input.ref} > {log} 2>&1

        # samtools faidx  -> .fai
        samtools faidx {input.ref} >> {log} 2>&1

        # GATK CreateSequenceDictionary -> .dict
        gatk CreateSequenceDictionary \
            -R {input.ref} \
            -O {output.dict} \
            >> {log} 2>&1

        touch {output.bwa_ok}
        """


# ----------------------------
# 1) fastp
# ----------------------------
rule fastp:
    input:
        refok=rules.ref_index.output,
        r1=lambda wc: _r1(wc.sample),
        r2=lambda wc: _r2(wc.sample) if PAIRED == 2 else None,
    output:
        r1c=f"{FASTP_DIR}/{{sample}}_R1.clean.fastq.gz",
        r2c=f"{FASTP_DIR}/{{sample}}_R2.clean.fastq.gz" if PAIRED == 2 else f"{FASTP_DIR}/{{sample}}_R2.clean.fastq.gz.dummy",
        json=f"{FASTP_DIR}/{{sample}}.fastp.json",
        html=f"{FASTP_DIR}/{{sample}}.fastp.html",
    threads: THREADS
    log:
        f"{LOG_DIR}/fastp/{{sample}}.log"
    benchmark:
        f"{BENCH_DIR}/fastp_{{sample}}.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {FASTP_DIR} {LOG_DIR}/fastp

        if [ "{PAIRED}" -eq 2 ]; then
          fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1c} -O {output.r2c} \
            -j {output.json} -h {output.html} \
            -w {threads} \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --trim_poly_x \
            > {log} 2>&1
          test -s {output.r2c}
        else
          fastp \
            -i {input.r1} \
            -o {output.r1c} \
            -j {output.json} -h {output.html} \
            -w {threads} \
            --trim_poly_g \
            --trim_poly_x \
            > {log} 2>&1
          echo "SE" > {output.r2c}
        fi

        test -s {output.r1c}
        test -s {output.json}
        test -s {output.html}
        """

# ----------------------------
# 2) MultiQC
# ----------------------------
rule multiqc:
    input:
        expand(f"{FASTP_DIR}/{{sample}}.fastp.json", sample=SAMPLES)
    output:
        f"{MULTIQC_DIR}/multiqc_report.html"
    log:
        f"{LOG_DIR}/multiqc.log"
    benchmark:
        f"{BENCH_DIR}/multiqc.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {MULTIQC_DIR}
        multiqc {FASTP_DIR} -o {MULTIQC_DIR} --force > {log} 2>&1
        test -s {output}
        """
        
# ----------------------------
# 3) BWA mem + sort + mapping stats (MAPPED, pre-dedup)
# ----------------------------
rule bwa_map:
    input:
        refok = rules.ref_index.output.bwa_ok,
        ref = REF_FASTA,
        r1  = f"{FASTP_DIR}/{{sample}}_R1.clean.fastq.gz",
        r2  = f"{FASTP_DIR}/{{sample}}_R2.clean.fastq.gz" if PAIRED == 2 else None,
    output:
        bam   = f"{BWA_DIR}/{{sample}}.sorted.bam",          # mapped (pre-dedup)
        bai   = f"{BWA_DIR}/{{sample}}.sorted.bam.bai",
        flag  = f"{BWA_DIR}/{{sample}}.flagstat.txt",
        idx   = f"{BWA_DIR}/{{sample}}.idxstats.txt",
        depth = f"{BWA_DIR}/{{sample}}.mapped.depth.tsv",    # mapped depth
    log:
        f"{LOG_DIR}/bwa/{{sample}}.log"
    benchmark:
        f"{BENCH_DIR}/bwa_{{sample}}.txt"
    threads: 4
    resources:
        tmpdir=lambda wc: f"{TMP_DIR}/{wc.sample}"

    run:
        import os
        from snakemake.shell import shell

        os.makedirs(BWA_DIR, exist_ok=True)
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        os.makedirs(resources.tmpdir, exist_ok=True)

        rg = (
            f"@RG\\tID:{wildcards.sample}"
            f"\\tSM:{wildcards.sample}"
            f"\\tPL:ILLUMINA"
            f"\\tLB:{wildcards.sample}"
            f"\\tPU:unit1"
        )

        # -----------------------------------
        # bwa mem → sorted BAM (mapped, pre-dedup)
        # -----------------------------------
        if PAIRED == 2:
            shell(
                r"""
                set -euo pipefail
                bash -c '
                  bwa mem -t {threads} -R "{rg}" \
                    {input.ref} {input.r1} {input.r2} \
                  | samtools sort -@ {threads} -m 1G \
                      -T {resources.tmpdir}/{wildcards.sample} \
                      -o {output.bam} -
                ' > {log} 2>&1
                """
            )
        else:
            shell(
                r"""
                set -euo pipefail
                bash -c '
                  bwa mem -t {threads} -R "{rg}" \
                    {input.ref} {input.r1} \
                  | samtools sort -@ {threads} -m 1G \
                      -T {resources.tmpdir}/{wildcards.sample} \
                      -o {output.bam} -
                ' > {log} 2>&1
                """
            )

        # -----------------------------------
        # index + basic mapping stats
        # -----------------------------------
        shell(f"samtools index {output.bam} >> {log} 2>&1")
        shell(f"samtools flagstat {output.bam} > {output.flag} 2>> {log}")
        shell(f"samtools idxstats {output.bam} > {output.idx} 2>> {log}")

        # -----------------------------------
        # depth (mapped, pre-dedup)
        # -----------------------------------
        shell(
            f"samtools depth -aa -d 0 "
            f"{output.bam} > {output.depth} 2>> {log}"
        )

        # sanity
        shell(f"test -s {output.bam}")
        shell(f"test -s {output.depth}")
            
# ----------------------------
# 4) Picard MarkDuplicates (PROCESSED BAM)
# ----------------------------
rule picard_markdup:
    input:
        bam = rules.bwa_map.output.bam   # mapped (pre-dedup)
    output:
        bam = f"{BAMQC_DIR}/{{sample}}.md.bam",          # processed BAM
        # 핵심 수정: .bam.bai -> .bai (Picard 생성 규칙 적용)
        bai = f"{BAMQC_DIR}/{{sample}}.md.bai",
        metrics = f"{BAMQC_DIR}/{{sample}}.markdup.metrics.txt"
    threads: 1
    resources:
        picard=1
    log:
        f"{LOG_DIR}/picard/{{sample}}.log"
    benchmark:
        f"{BENCH_DIR}/picard_markdup_{{sample}}.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {BAMQC_DIR} {LOG_DIR}/picard

        TMPDIR="{TMP_DIR}/picard_{wildcards.sample}_$$"
        mkdir -p "$TMPDIR"

        # -----------------------------------
        # MarkDuplicates 실행
        # -----------------------------------
          gatk --java-options "-Xmx8g -Djava.io.tmpdir=$TMPDIR" \
            MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT \
          > {log} 2>&1

        # -----------------------------------
        # sanity check (BAM 파일 손상 여부 확인)
        # -----------------------------------
        samtools quickcheck -v {output.bam} \
          >> {log} 2>&1

        rm -rf "$TMPDIR"
        
        # 파일이 실제로 생성되었고 크기가 0보다 큰지 최종 확인
        test -s {output.bam}
        test -s {output.bai}
        test -s {output.metrics}
        """

rule bam_sanity:
    input:
        # picard_markdup의 결과물을 입력으로 받음
        bam = f"{BAMQC_DIR}/{{sample}}.md.bam"
    output:
        ok  = f"{BAMQC_DIR}/{{sample}}.md.bam.ok"
    log:
        f"{LOG_DIR}/bam_sanity/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {LOG_DIR}/bam_sanity

        # samtools view로 실제 읽히는지 한 번 더 확인
        samtools view -c "{input.bam}" \
          > /dev/null 2> "{log}"

        touch "{output.ok}"
        """

# ----------------------------
# Calling을 위한 BAM/BAI 경로 결정 함수
# ----------------------------
def _bam_for_calling(sample):
    # Markdup 사용 여부에 따라 경로 반환
    if USE_MARKDUP == 1:
        return f"{BAMQC_DIR}/{sample}.md.bam"
    return f"{BWA_DIR}/{sample}.sorted.bam"

def _bai_for_calling(sample):
    if USE_MARKDUP == 1:
        # 핵심 수정: 위 규칙의 output과 일치하도록 .md.bai 로 변경
        return f"{BAMQC_DIR}/{sample}.md.bai"
    return f"{BWA_DIR}/{sample}.sorted.bam.bai"

# ----------------------------
# 5) GATK HaplotypeCaller → VCF.gz
# ----------------------------
rule gatk_call:
    input:
        refok = rules.ref_index.output.bwa_ok,
        ref   = REF_FASTA,
        bam   = lambda wc: _bam_for_calling(wc.sample),
        bai   = lambda wc: _bai_for_calling(wc.sample),
        mdok  = rules.picard_markdup.output.bam if USE_MARKDUP == 1 else rules.bwa_map.output.bam
    output:
        vcf = f"{GATK_DIR}/{{sample}}.gatk.vcf.gz",
        tbi = f"{GATK_DIR}/{{sample}}.gatk.vcf.gz.tbi",
    threads: THREADS
    log:
        f"{LOG_DIR}/gatk/{{sample}}.log"
    benchmark:
        f"{BENCH_DIR}/gatk_{{sample}}.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {GATK_DIR} {LOG_DIR}/gatk

        # -----------------------------------
        # GATK HaplotypeCaller (haploid)
        # -----------------------------------
        gatk --java-options "-Xmx8g" HaplotypeCaller \
          -R "{input.ref}" \
          -I "{input.bam}" \
          --sample-ploidy 1 \
          --min-base-quality-score 20 \
          --pcr-indel-model NONE \
          --native-pair-hmm-threads {threads} \
          -O "{output.vcf}" \
          2>&1 | tee "{log}"

        # -----------------------------------
        # tabix index
        # -----------------------------------
        tabix -f -p vcf "{output.vcf}" \
          2>&1 | tee -a "{log}"

        # -----------------------------------
        # sanity check
        # -----------------------------------
        test -s "{output.vcf}"
        test -s "{output.tbi}"
    """

# ----------------------------
# 6) Variant stats (bcftools stats)
# ----------------------------
rule bcftools_stats_gatk:
    input:
        vcf=rules.gatk_call.output.vcf,
        tbi=rules.gatk_call.output.tbi,
    output:
        stats=f"{SUM_GATK}/{{sample}}.bcftools.stats.txt"
    log:
        f"{LOG_DIR}/bcftools/gatk_{{sample}}.log"
    benchmark:
        f"{BENCH_DIR}/bcftools_gatk_{{sample}}.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {SUM_GATK} {LOG_DIR}/bcftools
        bcftools stats {input.vcf} > {output.stats} 2> {log}
        test -s {output.stats}
        """

# ----------------------------
# 7) CNV (window-based, PROCESSED BAM)
# ----------------------------
rule cnv_window:
    input:
        bam = lambda wc: _bam_for_calling(wc.sample),
        bai = lambda wc: _bai_for_calling(wc.sample),
    output:
        tsv = f"{CNV_WIN_DIR}/{{sample}}.cnv.window.tsv"
    log:
        f"{LOG_DIR}/cnv/{{sample}}.cnv_window.log"
    benchmark:
        f"{BENCH_DIR}/cnv_window_{{sample}}.txt"
    run:
        import pandas as pd
        import numpy as np
        import subprocess
        import os
        import tempfile

        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)

        # -----------------------------------
        # 1) depth from PROCESSED BAM
        # -----------------------------------
        with tempfile.NamedTemporaryFile(delete=False, suffix=".depth") as tmp:
            depth_file = tmp.name

        subprocess.run(
            [
                "samtools", "depth", "-aa", "-d", "0",
                input.bam
            ],
            stdout=open(depth_file, "w"),
            stderr=open(log[0], "w"),
            check=True
        )

        df = pd.read_csv(
            depth_file,
            sep="\t",
            header=None,
            names=["chrom", "pos", "depth"]
        )

        os.unlink(depth_file)

        if df.empty:
            pd.DataFrame(
                columns=["chrom","start","end","mean_depth","log2_ratio"]
            ).to_csv(output.tsv, sep="\t", index=False)
            return

        # -----------------------------------
        # 2) window aggregation
        # -----------------------------------
        df["win_start"] = ((df["pos"] - 1) // CNV_WINDOW) * CNV_WINDOW + 1
        df["win_end"]   = df["win_start"] + CNV_WINDOW - 1

        g = (
            df
            .groupby(["chrom","win_start","win_end"], as_index=False)["depth"]
            .mean()
            .rename(columns={
                "win_start":"start",
                "win_end":"end",
                "depth":"mean_depth"
            })
        )

        # -----------------------------------
        # 3) normalize by genome-wide median
        # -----------------------------------
        genome_median = float(np.median(df["depth"].values))
        genome_median = genome_median if genome_median > 0 else 1.0

        g["log2_ratio"] = np.log2((g["mean_depth"] + 1e-6) / genome_median)

        # -----------------------------------
        # 4) write output
        # -----------------------------------
        g.to_csv(output.tsv, sep="\t", index=False)

        with open(log[0], "a") as fh:
            fh.write(
                f"CNV_WINDOW={CNV_WINDOW}\n"
                f"processed_bam={input.bam}\n"
                f"genome_median_depth={genome_median}\n"
            )

# ----------------------------
# 8) CNV (gene-based, PROCESSED BAM)
# ----------------------------
rule cnv_gene:
    input:
        bam = lambda wc: _bam_for_calling(wc.sample),
        bai = lambda wc: _bai_for_calling(wc.sample),
    output:
        tsv = f"{CNV_GENE_DIR}/{{sample}}.cnv.gene.tsv"
    log:
        f"{LOG_DIR}/cnv/{{sample}}.cnv_gene.log"
    benchmark:
        f"{BENCH_DIR}/cnv_gene_{{sample}}.txt"
    run:
        import pandas as pd
        import numpy as np
        import subprocess
        import tempfile
        import os
        import re

        if not GFF_FILE:
            raise ValueError("[FATAL] cnv_gene requires config: gff=/path/to/PlasmoDB.gff")

        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)

        # -----------------------------------
        # 1) depth from PROCESSED BAM
        # -----------------------------------
        with tempfile.NamedTemporaryFile(delete=False, suffix=".depth") as tmp:
            depth_file = tmp.name

        subprocess.run(
            [
                "samtools", "depth", "-aa", "-d", "0",
                input.bam
            ],
            stdout=open(depth_file, "w"),
            stderr=open(log[0], "w"),
            check=True
        )

        df_depth = pd.read_csv(
            depth_file,
            sep="\t",
            header=None,
            names=["chrom", "pos", "depth"]
        )
        os.unlink(depth_file)

        if df_depth.empty:
            pd.DataFrame(
                columns=[
                    "sample","gene_id","chrom","start","end",
                    "strand","mean_depth","log2_ratio"
                ]
            ).to_csv(output.tsv, sep="\t", index=False)
            return

        # genome-wide median depth
        genome_median = float(np.median(df_depth["depth"].values))
        genome_median = genome_median if genome_median > 0 else 1.0

        # -----------------------------------
        # 2) build prefix sums per chromosome
        # -----------------------------------
        chrom_ps = {}
        for chrom, sub in df_depth.groupby("chrom"):
            maxpos = int(sub["pos"].max())
            arr = np.zeros(maxpos + 1, dtype=np.float32)
            arr[sub["pos"].values.astype(int)] = sub["depth"].values.astype(np.float32)
            chrom_ps[chrom] = np.cumsum(arr)

        # -----------------------------------
        # 3) parse GFF (gene features)
        # -----------------------------------
        genes = []
        with open(GFF_FILE) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                chrom, source, feature, start, end, score, strand, phase, attrs = parts
                if feature not in ["gene", "protein_coding_gene", "ncRNA_gene", "pseudogene"]:
                    continue
                m = re.search(r"(?:^|;)ID=([^;]+)", attrs)
                if not m:
                    m = re.search(r"(?:^|;)gene_id=([^;]+)", attrs)
                if not m:
                    continue

                genes.append({
                    "gene_id": m.group(1),
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand
                })

        # -----------------------------------
        # 4) gene-level mean depth + log2 ratio
        # -----------------------------------
        rows = []
        for g in genes:
            chrom = g["chrom"]
            if chrom not in chrom_ps:
                continue

            ps = chrom_ps[chrom]
            s = max(1, g["start"])
            e = min(len(ps) - 1, g["end"])
            if e < s:
                continue

            total = ps[e] - ps[s - 1]
            mean_depth = total / float(e - s + 1)
            log2_ratio = float(np.log2((mean_depth + 1e-6) / genome_median))

            rows.append({
                "sample": wildcards.sample,
                "gene_id": g["gene_id"],
                "chrom": chrom,
                "start": s,
                "end": e,
                "strand": g["strand"],
                "mean_depth": mean_depth,
                "log2_ratio": log2_ratio
            })

        df_out = pd.DataFrame(
            rows,
            columns=[
                "sample","gene_id","chrom","start","end",
                "strand","mean_depth","log2_ratio"
            ]
        )

        df_out.to_csv(output.tsv, sep="\t", index=False)

        with open(log[0], "a") as fh:
            fh.write(
                f"GFF={GFF_FILE}\n"
                f"processed_bam={input.bam}\n"
                f"genome_median_depth={genome_median}\n"
            )


# ----------------------------
# 9) Mixed infection metrics (AF spectrum + heterozygosity proxy)
# ----------------------------
rule af_table:
    input:
        vcf=rules.gatk_call.output.vcf,
        tbi=rules.gatk_call.output.tbi
    output:
        tsv=f"{MIX_DIR}/{{sample}}.af.tsv"
    log:
        f"{LOG_DIR}/mixed/{{sample}}.af.log"
    benchmark:
        f"{BENCH_DIR}/mixed_af_{{sample}}.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {MIX_DIR} {LOG_DIR}/mixed
        # Extract AD and DP if present. AD is typically "ref,alt" for haploid too.
        # Output columns: CHROM POS REF ALT AD DP
        bcftools query \
          -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\t%DP\n' \
          {input.vcf} > {output.tsv} 2> {log}
        test -s {output.tsv}
        """

rule mixed_metrics:
    input:
        af=rules.af_table.output.tsv
    output:
        tsv=f"{MIX_DIR}/{{sample}}.mixed.metrics.tsv"
    log:
        f"{LOG_DIR}/mixed/{{sample}}.mixed_metrics.log"
    benchmark:
        f"{BENCH_DIR}/mixed_metrics_{{sample}}.txt"
    run:
        import pandas as pd, numpy as np

        rows = []
        with open(input.af) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                chrom, pos, ref, alt, ad, dp = parts[:6]
                # AD might be like "12,3" or "."; DP might be "."
                if ad == "." or "," not in ad:
                    continue
                try:
                    refc, altc = ad.split(",", 1)
                    refc = float(refc); altc = float(altc)
                    dpv = refc + altc
                    if dpv <= 0:
                        continue
                    af = altc / dpv
                except Exception:
                    continue
                rows.append(af)

        if not rows:
            # empty -> write minimal
            out = pd.DataFrame([{
                "sample": wildcards.sample,
                "n_sites_with_AD": 0,
                "mean_AF": None,
                "frac_intermediate_AF_0.2_0.8": None,
                "heterozygosity_proxy_mean_2p1p": None
            }])
            out.to_csv(output.tsv, sep="\t", index=False)
            open(log[0], "w").write("No AD sites found\n")
            return

        afs = np.array(rows, dtype=float)
        het = 2.0 * afs * (1.0 - afs)  # proxy (max at 0.5)
        frac_mid = float(((afs >= 0.2) & (afs <= 0.8)).mean())
        out = pd.DataFrame([{
            "sample": wildcards.sample,
            "n_sites_with_AD": int(len(afs)),
            "mean_AF": float(afs.mean()),
            "frac_intermediate_AF_0.2_0.8": frac_mid,
            "heterozygosity_proxy_mean_2p1p": float(het.mean())
        }])
        out.to_csv(output.tsv, sep="\t", index=False)
        open(log[0], "w").write("Computed AF-based mixed infection metrics\n")

# ----------------------------
# 10) Recombination exploratory (PROCESSED BAM)
#     - subtelomeric soft-clipped reads
#     - discordant read pairs (PE only)
# ----------------------------
rule recomb_explore:
    input:
        bam = lambda wc: _bam_for_calling(wc.sample),   # processed BAM
        bai = lambda wc: _bai_for_calling(wc.sample),
        fai = REF_FASTA + ".fai"
    output:
        tsv = f"{RECOMB_DIR}/{{sample}}.recomb_explore.tsv"
    log:
        f"{LOG_DIR}/recomb/{{sample}}.recomb.log"
    benchmark:
        f"{BENCH_DIR}/recomb_explore_{{sample}}.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {RECOMB_DIR} {LOG_DIR}/recomb

        # ------------------------------------------------------------
        # Build subtelomeric BED (from .fai)
        # ------------------------------------------------------------
        BED_TMP="$(mktemp -p /tmp {wildcards.sample}.subtel.XXXXXX.bed)"
        awk -v S={SUBTEL_SIZE} 'BEGIN{{OFS="\t"}} {{
            L=$2;
            if (L > S) {{
                print $1, 0, S;
                print $1, L-S, L;
            }} else {{
                print $1, 0, L;
            }}
        }}' {input.fai} > "$BED_TMP"

        # ------------------------------------------------------------
        # Soft-clipped reads (processed BAM)
        # ------------------------------------------------------------
        SOFTCLIP=$(
          while read -r c s e; do
            samtools view {input.bam} "$c:$((s+1))-$e" 2>> {log} \
              | awk '$6 ~ /S/ {{n++}} END{{print n+0}}'
          done < "$BED_TMP" | awk '{{sum+=$1}} END{{print sum+0}}'
        )

        # ------------------------------------------------------------
        # Discordant pairs (processed BAM, PE only)
        # ------------------------------------------------------------
        if [ "{PAIRED}" -eq 2 ]; then
          DISCORD=$(
            while read -r c s e; do
              samtools view \
                -f 1 -F 2 -F 4 -F 8 {input.bam} "$c:$((s+1))-$e" 2>> {log} \
                | wc -l
            done < "$BED_TMP" | awk '{{sum+=$1}} END{{print sum+0}}'
          )
        else
          DISCORD="NA"
        fi

        # ------------------------------------------------------------
        # Write output (atomic)
        # ------------------------------------------------------------
        OUTTMP="$(mktemp -p /tmp {wildcards.sample}.recomb.XXXXXX.tsv)"
        echo -e "sample\tSUBTEL_SIZE\tsoftclipped_reads\tdiscordant_pairs" > "$OUTTMP"
        echo -e "{wildcards.sample}\t{SUBTEL_SIZE}\t$SOFTCLIP\t$DISCORD" >> "$OUTTMP"
        mv -f "$OUTTMP" "{output.tsv}"

        rm -f "$BED_TMP"
        test -s "{output.tsv}"
        """

# ----------------------------
# 11) Drug resistance gene summary (requires GFF)
#     - extract gene coordinates from GFF (by gene ID)
#     - bcftools view -r region to get variants
# ----------------------------
rule drug_variants:
    input:
        vcf=rules.gatk_call.output.vcf,
        tbi=rules.gatk_call.output.tbi
    output:
        tsv=f"{DRUG_DIR}/{{sample}}.drug_variants.tsv"
    log:
        f"{LOG_DIR}/drug/{{sample}}.drug.log"
    benchmark:
        f"{BENCH_DIR}/drug_{{sample}}.txt"
    run:
        import pandas as pd
        import subprocess
        import re

        if not GFF_FILE:
            raise ValueError("[FATAL] Drug summary requires config: gff=/path/to/PlasmoDB.gff")

        wanted = set(DRUG_GENES.values())
        coords = {}  # gene_id -> (chrom, start, end)

        with open(GFF_FILE) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue

                chrom, source, feature, start, end, score, strand, phase, attrs = parts

                # 🔧 핵심 수정 1: PlasmoDB gene feature 허용
                if feature not in ("gene", "protein_coding_gene", "pseudogene", "ncRNA_gene"):
                    continue

                m = re.search(r"(?:^|;)ID=([^;]+)", attrs)
                if not m:
                    continue

                gid = m.group(1)

                # 🔧 핵심 수정 2: transcript suffix 제거 (PF3D7_xxx.1 → PF3D7_xxx)
                gid = gid.split(".")[0]

                if gid in wanted:
                    coords[gid] = (chrom, int(start), int(end))

        rows = []
        for label, gid in DRUG_GENES.items():
            if gid not in coords:
                rows.append({
                    "sample": wildcards.sample,
                    "gene_label": label,
                    "gene_id": gid,
                    "region": None,
                    "n_variants": 0,
                    "note": "gene_not_found_in_gff"
                })
                continue

            chrom, s, e = coords[gid]
            region = f"{chrom}:{s}-{e}"

            cmd = [
                "bcftools", "query",
                "-r", region,
                "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\t%QUAL\n",
                input.vcf
            ]

            res = subprocess.run(cmd, capture_output=True, text=True)
            out = res.stdout.strip().splitlines() if res.stdout else []

            if not out:
                rows.append({
                    "sample": wildcards.sample,
                    "gene_label": label,
                    "gene_id": gid,
                    "region": region,
                    "n_variants": 0,
                    "variants": ""
                })
            else:
                rows.append({
                    "sample": wildcards.sample,
                    "gene_label": label,
                    "gene_id": gid,
                    "region": region,
                    "n_variants": len(out),
                    "variants": " | ".join(out[:500])
                })

        df = pd.DataFrame(rows)
        df.to_csv(output.tsv, sep="\t", index=False)

        with open(log[0], "w") as f:
            f.write(f"GFF={GFF_FILE}\nDRUG_GENES={DRUG_GENES}\n")


# -----------------------------------
# Per-sample SNP table (for long → wide pivot)
# -----------------------------------
rule snp_table_per_sample:
    input:
        vcf = rules.gatk_call.output.vcf,
        tbi = rules.gatk_call.output.tbi
    output:
        tsv = f"{SUM_DIR}/snp_table/{{sample}}.snps.tsv"
    log:
        f"{LOG_DIR}/snp/{{sample}}.snps.log"
    benchmark:
        f"{BENCH_DIR}/snp_table_{{sample}}.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {SUM_DIR}/snp_table {LOG_DIR}/snp

        # SNP only (REF, ALT length == 1)
        # Extract ALT / DP / AF explicitly (PF haploid-friendly)
        bcftools query \
          -i 'strlen(REF)=1 && strlen(ALT)=1' \
          -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t[%AF]\n' \
          {input.vcf} > {output.tsv} 2> {log}

        # keep header if empty (safe for downstream merge)
        if [ ! -s {output.tsv} ]; then
          echo -e "chrom\tpos\tref\talt\tdp\taf" > {output.tsv}
        fi
        """
       
# -----------------------------------
# qc_report  (mapped vs processed)
# -----------------------------------
rule qc_report:
    input:
        mapped_bam = f"{BWA_DIR}/{{sample}}.sorted.bam",
        mapped_bai = f"{BWA_DIR}/{{sample}}.sorted.bam.bai",
        proc_bam   = lambda wc: _bam_for_calling(wc.sample),
        proc_bai   = lambda wc: _bai_for_calling(wc.sample),
        markdup_metrics = f"{BAMQC_DIR}/{{sample}}.markdup.metrics.txt",
    output:
        main = f"{QC_DIR}/{{sample}}.qc_report.tsv",
        chr  = f"{QC_DIR}/{{sample}}.qc_report.chr.tsv"
    log:
        f"{LOG_DIR}/qc/{{sample}}.qc_report.log"
    benchmark:
        f"{BENCH_DIR}/qc_report_{{sample}}.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {QC_DIR} {LOG_DIR}/qc {TMP_DIR}

        MBAM="{input.mapped_bam}"
        PBAM="{input.proc_bam}"
        TMP="{TMP_DIR}/{wildcards.sample}"
        MAIN_TMP="$TMP.main.tsv"
        CHR_TMP="$TMP.chr.tsv"

        ########################################
        # 1) samtools stats (processed BAM 기준)
        ########################################
        samtools stats "$PBAM" > "$TMP.stats" 2>> {log}

        read_mean=$(awk -F'\t' '$1=="SN" && $2~"average length" {{print $3; exit}}' "$TMP.stats")
        read_min=$(awk -F'\t' '$1=="SN" && $2~"minimum length" {{print $3; exit}}' "$TMP.stats")
        read_max=$(awk -F'\t' '$1=="SN" && $2~"maximum length" {{print $3; exit}}' "$TMP.stats")

        mean_mapq=$(awk -F'\t' '$1=="SN" && $2~"average quality" {{print $3; exit}}' "$TMP.stats")
        error_rate=$(awk -F'\t' '$1=="SN" && $2~"error rate" {{print $3; exit}}' "$TMP.stats")

        ########################################
        # 2) read counts (mapped vs processed)
        ########################################
        raw_reads=$(samtools view -c -F 256 -F 2048 "$MBAM")
        mapped_reads=$(samtools view -c -F 4 -F 256 -F 2048 "$MBAM")
        processed_reads=$(samtools view -c -F 4 -F 256 -F 2048 "$PBAM")



        mapped_pct=$(awk -v m=$mapped_reads -v t=$raw_reads \
          'BEGIN{{printf "%.2f", (t>0)?(m/t*100):0}}')

        dup_reads=$(awk '
          BEGIN{{FS="\t"}}
          $1=="LIBRARY"{{
            getline;
            print $6 + $7;
            exit
          }}
        ' {input.markdup_metrics})

        dup_rate=$(awk '
          BEGIN{{FS="\t"}}
          $1=="LIBRARY"{{
            getline;
            printf "%.6f", $9;
            exit
          }}
        ' {input.markdup_metrics})

        ########################################
        # 3) depth: mapped (pre-dedup)
        ########################################
        samtools depth -aa -d 0 "$MBAM" > "$TMP.mapped.depth"

        mapped_stats=$(awk '
          BEGIN{{sum=0;sumsq=0;n=0;c10=0;c20=0;c30=0}}
          {{
            d=$3;
            sum+=d; sumsq+=d*d; n++;
            if(d>=10) c10++;
            if(d>=20) c20++;
            if(d>=30) c30++;
          }}
          END{{
            mean=sum/n;
            sd=sqrt((sumsq/n)-(mean*mean));
            printf "%.4f\t%.4f\t%.2f\t%.2f\t%.2f\n",
              mean, sd, c10/n*100, c20/n*100, c30/n*100
          }}
        ' "$TMP.mapped.depth")

        mapped_mean=$(echo "$mapped_stats" | awk '{{print $1}}')
        mapped_sd=$(echo "$mapped_stats"   | awk '{{print $2}}')
        mapped_10x=$(echo "$mapped_stats"  | awk '{{print $3}}')
        mapped_20x=$(echo "$mapped_stats"  | awk '{{print $4}}')
        mapped_30x=$(echo "$mapped_stats"  | awk '{{print $5}}')

        ########################################
        # 4) depth: processed (post-dedup)
        ########################################
        samtools depth -aa -d 0 "$PBAM" > "$TMP.proc.depth"

        proc_stats=$(awk '
          BEGIN{{sum=0;sumsq=0;n=0;c10=0;c20=0;c30=0}}
          {{
            d=$3;
            sum+=d; sumsq+=d*d; n++;
            if(d>=10) c10++;
            if(d>=20) c20++;
            if(d>=30) c30++;
          }}
          END{{
            mean=sum/n;
            sd=sqrt((sumsq/n)-(mean*mean));
            printf "%.4f\t%.4f\t%.2f\t%.2f\t%.2f\n",
              mean, sd, c10/n*100, c20/n*100, c30/n*100
          }}
        ' "$TMP.proc.depth")

        proc_mean=$(echo "$proc_stats" | awk '{{print $1}}')
        proc_sd=$(echo "$proc_stats"   | awk '{{print $2}}')
        proc_10x=$(echo "$proc_stats"  | awk '{{print $3}}')
        proc_20x=$(echo "$proc_stats"  | awk '{{print $4}}')
        proc_30x=$(echo "$proc_stats"  | awk '{{print $5}}')

        ########################################
        # 5) per-chrom depth (processed 기준)
        ########################################
        awk '
          BEGIN{{OFS="\t"; print "chrom","mean_depth","sd_depth"}}
          {{
            c=$1; d=$3;
            sum[c]+=d; sumsq[c]+=d*d; n[c]++
          }}
          END{{
            for(c in n){{
              m=sum[c]/n[c];
              sd=sqrt((sumsq[c]/n[c])-(m*m));
              printf "%s\t%.4f\t%.4f\n", c, m, sd
            }}
          }}
        ' "$TMP.proc.depth" | sort -k1,1 > "$CHR_TMP"

        ########################################
        # 6) main QC table
        ########################################
        echo -e \
"sample\traw_reads\tmapped_reads\tprocessed_reads\tmapped_pct\tmapped_mean_depth\tmapped_sd_depth\tmapped_10X\tmapped_20X\tmapped_30X\tprocessed_mean_depth\tprocessed_sd_depth\tprocessed_10X\tprocessed_20X\tprocessed_30X\tdup_reads\tduplication_rate\tread_min\tread_max\tread_mean\tmean_mapq\terror_rate" \
> "$MAIN_TMP"

        echo -e \
"{wildcards.sample}\t$raw_reads\t$mapped_reads\t$processed_reads\t$mapped_pct\t$mapped_mean\t$mapped_sd\t$mapped_10x\t$mapped_20x\t$mapped_30x\t$proc_mean\t$proc_sd\t$proc_10x\t$proc_20x\t$proc_30x\t$dup_reads\t$dup_rate\t$read_min\t$read_max\t$read_mean\t$mean_mapq\t$error_rate" \
>> "$MAIN_TMP"

        mv -f "$MAIN_TMP" "{output.main}"
        mv -f "$CHR_TMP"  "{output.chr}"
        rm -f "$TMP."*
        """

# -----------------------------------
# Merge per-sample SNP tables (LONG format)
# -----------------------------------
rule snp_table_all:
    input:
        snps = lambda wc: expand(
            f"{SUM_DIR}/snp_table/{{sample}}.snps.tsv",
            sample=SAMPLES
        )
    output:
        tsv = f"{SUM_DIR}/all_samples.snps.long.tsv"
    log:
        f"{LOG_DIR}/snp/snp_merge.log"
    benchmark:
        f"{BENCH_DIR}/snp_table_all_long.txt"
    run:
        import pandas as pd
        import os

        dfs = []

        for f in input.snps:
            sample = os.path.basename(f).replace(".snps.tsv", "")
            if os.path.getsize(f) == 0:
                continue

            df = pd.read_csv(
                f,
                sep="\t",
                names=["chrom","pos","ref","alt","dp","af"]
            )

            # empty / header-only safety
            if df.empty or df.iloc[0]["chrom"] == "chrom":
                continue

            df.insert(0, "sample", sample)
            dfs.append(df)

        if dfs:
            out = (
                pd.concat(dfs, ignore_index=True)
                .sort_values(["chrom","pos","sample"])
                .reset_index(drop=True)
            )
        else:
            out = pd.DataFrame(
                columns=["sample","chrom","pos","ref","alt","dp","af"]
            )

        out.to_csv(output.tsv, sep="\t", index=False)

# -----------------------------------
# SNP table (WIDE format for reporting)
# chrom | pos | ref | sampleA_alt | sampleA_dp | sampleA_af | sampleB_...
# -----------------------------------
rule snp_table_wide_full:
    input:
        long = f"{SUM_DIR}/all_samples.snps.long.tsv"
    output:
        tsv  = f"{SUM_DIR}/all_samples.snps.wide.tsv"
    log:
        f"{LOG_DIR}/snp/snp_wide.log"
    benchmark:
        f"{BENCH_DIR}/snp_table_wide_full.txt"
    run:
        import pandas as pd

        df = pd.read_csv(input.long, sep="\t")

        if df.empty:
            pd.DataFrame(
                columns=["chrom","pos","ref"]
            ).to_csv(output.tsv, sep="\t", index=False)
            return

        samples = sorted(df["sample"].unique())

        # 기본 좌표
        base = (
            df[["chrom","pos","ref"]]
            .drop_duplicates()
            .sort_values(["chrom","pos"])
            .reset_index(drop=True)
        )

        out = base.copy()

        # 🔥 sample 단위로 alt / dp / af를 연속으로 붙임
        for s in samples:
            sub = df[df["sample"] == s][["chrom","pos","ref","alt","dp","af"]]

            out = (
                out
                .merge(
                    sub.rename(columns={
                        "alt": f"{s}_alt",
                        "dp":  f"{s}_dp",
                        "af":  f"{s}_af",
                    }),
                    on=["chrom","pos","ref"],
                    how="left"
                )
            )

        out.to_csv(output.tsv, sep="\t", index=False)



# ----------------------------
# Final QC XLSX (mapped vs processed)
# ----------------------------
rule final_qc_xlsx:
    input:
        multiqc = f"{MULTIQC_DIR}/multiqc_report.html",
        fastp_jsons = lambda wc: expand(
            f"{FASTP_DIR}/{{sample}}.fastp.json", sample=SAMPLES
        ),
        qc_reports = lambda wc: expand(
            f"{QC_DIR}/{{sample}}.qc_report.tsv", sample=SAMPLES
        ),
    output:
        xlsx = FINAL_QC_XLSX
    log:
        f"{LOG_DIR}/final_qc_xlsx.log"
    benchmark:
        f"{BENCH_DIR}/final_qc_xlsx.txt"
    run:
        import pandas as pd, json, os

        # ==================================================
        # 1) fastp summary (RAW → CLEAN, read-level)
        # ==================================================
        fastp_rows = []
        for jp in input.fastp_jsons:
            with open(jp) as fh:
                d = json.load(fh)

            sample = os.path.basename(jp).replace(".fastp.json", "")
            bef = d["summary"]["before_filtering"]
            aft = d["summary"]["after_filtering"]

            fastp_rows.append({
                "sample": sample,
                "raw_reads_fastp": bef["total_reads"],
                "clean_reads_fastp": aft["total_reads"],
                "q30_rate": aft["q30_rate"],
                "gc_content": aft["gc_content"],
            })

        df_fastp = (
            pd.DataFrame(fastp_rows)
            .sort_values("sample")
            .reset_index(drop=True)
        )

        # ==================================================
        # 2) QC report (AUTHORITATIVE: mapped + processed)
        # ==================================================
        qc_dfs = [pd.read_csv(p, sep="\t") for p in input.qc_reports]
        df_qc = (
            pd.concat(qc_dfs, ignore_index=True)
            .sort_values("sample")
            .reset_index(drop=True)
        )

        # ==================================================
        # 3) Sample summary table (센터 스타일)
        # ==================================================
        df_summary = (
            df_qc[[
                "sample",

                # read counts
                "raw_reads",
                "mapped_reads",
                "processed_reads",
                "mapped_pct",

                # mapped (pre-dedup)
                "mapped_mean_depth",
                "mapped_sd_depth",
                "mapped_10X",
                "mapped_20X",
                "mapped_30X",

                # processed (post-dedup)
                "processed_mean_depth",
                "processed_sd_depth",
                "processed_10X",
                "processed_20X",
                "processed_30X",

                # duplicates / quality
                "dup_reads",
                "duplication_rate",
                "mean_mapq",
                "error_rate",
            ]]
            .merge(df_fastp, on="sample", how="left")
        )

        # ==================================================
        # 4) Write Excel
        # ==================================================
        os.makedirs(os.path.dirname(output.xlsx), exist_ok=True)
        with pd.ExcelWriter(output.xlsx, engine="openpyxl") as xw:

            # 핵심 요약
            df_summary.to_excel(
                xw, sheet_name="sample_summary", index=False
            )

            # QC 전체 컬럼
            df_qc.to_excel(
                xw, sheet_name="qc_report_full", index=False
            )

            # fastp (read-level)
            df_fastp.to_excel(
                xw, sheet_name="fastp", index=False
            )

            # meta
            pd.DataFrame({
                "run_tag": [RUN_TAG],
                "paired": [PAIRED],
                "use_markdup": [USE_MARKDUP],
                "processed_definition": ["Picard MarkDuplicates BAM"],
                "multiqc_html": [input.multiqc],
            }).to_excel(
                xw, sheet_name="meta", index=False
            )
            
# ----------------------------
# Final Variant XLSX
#  - variant summary
#  - CNV
#  - mixed infection
#  - recombination
#  - drug resistance
#  - SNPs (WIDE = main, LONG = optional)
# ----------------------------
rule final_variant_xlsx:
    input:
        gatk_stats = expand(
            f"{SUM_GATK}/{{sample}}.bcftools.stats.txt",
            sample=SAMPLES
        ),
        cnv_win = expand(
            f"{CNV_WIN_DIR}/{{sample}}.cnv.window.tsv",
            sample=SAMPLES
        ),
        cnv_gene = expand(
            f"{CNV_GENE_DIR}/{{sample}}.cnv.gene.tsv",
            sample=SAMPLES
        ) if GFF_FILE else [],
        mixed = expand(
            f"{MIX_DIR}/{{sample}}.mixed.metrics.tsv",
            sample=SAMPLES
        ),
        recomb = expand(
            f"{RECOMB_DIR}/{{sample}}.recomb_explore.tsv",
            sample=SAMPLES
        ),
        drug = expand(
            f"{DRUG_DIR}/{{sample}}.drug_variants.tsv",
            sample=SAMPLES
        ) if GFF_FILE else [],
        snps_long = f"{SUM_DIR}/all_samples.snps.long.tsv",
        snps_wide = f"{SUM_DIR}/all_samples.snps.wide.tsv"
    output:
        xlsx = FINAL_VAR_XLSX
    log:
        f"{LOG_DIR}/final_variant_xlsx.log"
    benchmark:
        f"{BENCH_DIR}/final_variant_xlsx.txt"
    run:
        import pandas as pd
        import os

        # -----------------------------------
        # 1) bcftools stats (summary)
        # -----------------------------------
        def parse_stats(path):
            sample = os.path.basename(path).replace(".bcftools.stats.txt","")
            d = {"sample": sample, "caller": "GATK"}

            with open(path) as f:
                for line in f:
                    if line.startswith("SN\t"):
                        p = line.strip().split("\t")
                        if len(p) >= 4:
                            d[p[2].rstrip(":")] = p[3]

                    elif line.startswith("TSTV\t"):
                        # TSTV  <multi> <Ts> <Tv> <Ts/Tv> <Ts_filt> <Tv_filt> <Ts/Tv_filt>
                        p = line.strip().split("\t")
                        if len(p) >= 8:
                            d["tstv_multiallelic"]   = int(p[1])
                            d["tstv_ts"]             = int(p[2])
                            d["tstv_tv"]             = int(p[3])
                            d["tstv_ratio"]          = float(p[4])
                            d["tstv_ts_filtered"]    = int(p[5])
                            d["tstv_tv_filtered"]    = int(p[6])
                            d["tstv_ratio_filtered"] = float(p[7])

            return d



        df_stats = (
            pd.DataFrame([parse_stats(p) for p in input.gatk_stats])
            .sort_values("sample")
        )

        keep = [
            "sample","caller",
            "number of records",
            "number of SNPs",
            "number of indels",
            "number of multiallelic sites",

            "tstv_multiallelic",
            "tstv_ts",
            "tstv_tv",
            "tstv_ratio",
            "tstv_ts_filtered",
            "tstv_tv_filtered",
            "tstv_ratio_filtered",
        ]
        df_stats_sum = df_stats[[c for c in keep if c in df_stats.columns]]

        # -----------------------------------
        # 2) Mixed infection
        # -----------------------------------
        df_mixed = (
            pd.concat([pd.read_csv(p, sep="\t") for p in input.mixed])
            .sort_values("sample")
        )

        # -----------------------------------
        # 3) Recombination explore
        # -----------------------------------
        df_recomb = (
            pd.concat([pd.read_csv(p, sep="\t") for p in input.recomb])
            .sort_values("sample")
        )

        # -----------------------------------
        # 4) CNV window summary
        # -----------------------------------
        cnv_rows = []
        for p in input.cnv_win:
            sample = os.path.basename(p).replace(".cnv.window.tsv","")
            d = pd.read_csv(p, sep="\t")
            cnv_rows.append({
                "sample": sample,
                "n_windows": len(d),
                "n_extreme_abs_log2_gt1": int((d["log2_ratio"].abs() > 1).sum()) if not d.empty else 0
            })
        df_cnv_win = pd.DataFrame(cnv_rows).sort_values("sample")

        # -----------------------------------
        # 5) CNV gene summary (optional)
        # -----------------------------------
        df_cnv_gene = None
        if GFF_FILE and input.cnv_gene:
            gene_rows = []
            for p in input.cnv_gene:
                sample = os.path.basename(p).replace(".cnv.gene.tsv","")
                d = pd.read_csv(p, sep="\t")
                gene_rows.append({
                    "sample": sample,
                    "n_genes": len(d),
                    "n_extreme_abs_log2_gt1": int((d["log2_ratio"].abs() > 1).sum()) if not d.empty else 0
                })
            df_cnv_gene = pd.DataFrame(gene_rows).sort_values("sample")

        # -----------------------------------
        # 6) Drug resistance (optional)
        # -----------------------------------
        df_drug = None
        if GFF_FILE and input.drug:
            df_drug = pd.concat(
                [pd.read_csv(p, sep="\t") for p in input.drug],
                ignore_index=True
            )

        # -----------------------------------
        # 7) SNP tables
        # -----------------------------------
        df_snps_wide = pd.read_csv(input.snps_wide, sep="\t")
        df_snps_long = pd.read_csv(input.snps_long, sep="\t")

        # -----------------------------------
        # 8) Write Excel
        # -----------------------------------
        os.makedirs(os.path.dirname(output.xlsx), exist_ok=True)
        with pd.ExcelWriter(output.xlsx, engine="openpyxl") as xw:

            df_stats_sum.to_excel(
                xw, sheet_name="variant_summary", index=False
            )
            df_stats.to_excel(
                xw, sheet_name="gatk_stats_full", index=False
            )

            df_mixed.to_excel(
                xw, sheet_name="mixed_metrics", index=False
            )
            df_recomb.to_excel(
                xw, sheet_name="recomb_explore", index=False
            )

            df_cnv_win.to_excel(
                xw, sheet_name="cnv_window_summary", index=False
            )
            if df_cnv_gene is not None:
                df_cnv_gene.to_excel(
                    xw, sheet_name="cnv_gene_summary", index=False
                )

            if df_drug is not None:
                df_drug.to_excel(
                    xw, sheet_name="drug_resistance", index=False
                )

            # ⭐ SNPs (main)
            df_snps_wide.to_excel(
                xw, sheet_name="snps_wide", index=False
            )

            # SNPs (debug / reference)
            # df_snps_long.to_excel(
            #     xw, sheet_name="snps_long", index=False
            # )

            # meta
            pd.DataFrame({
                "run_tag":[RUN_TAG],
                "processed_definition":["Picard MarkDuplicates BAM"],
                "snp_format_main":["wide"],
                "snp_format_aux":["long"]
            }).to_excel(
                xw, sheet_name="meta", index=False
            )
