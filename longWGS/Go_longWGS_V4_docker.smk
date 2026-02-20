#############################################
# Go_autocycler.smk
# 20260213
# Heekuk Park
#
# Pipeline:
#   QC → autocycler → medaka → quast → checkm2 → coverage → bakta → bandage
#
# OUTDIR/
# └── 1_QC/             (porechop_abi + NanoFilt output: *.clean.fastq.gz)
# └── 2_quast/          (QUAST reports; from medaka fasta)
# └── 3_autocycler/     ({sample}/autocycler_out/consensus_assembly.{fasta,gfa})
# └── 4_medaka/         ({sample}_final_assembly.fasta)
# └── 5_checkm2/        ({sample}/quality_report.tsv + DONE.txt + raw/)
# └── 6_coverage/       ({sample}.sorted.bam + {sample}_contig_mean_depth.tsv + DONE.txt)
# └── 7_bakta/          ({sample}/DONE.txt + outputs)
#
# DB structure (single root):
#   DB_ROOT/
#     ├── bakta_DB/
#     └── CheckM2_database/
#
# Example run:
# dir="test"
# output_dir="3_assembled_test"
# DB="/media/uhlemannlab/Nook01/DB"
#
# snakemake --snakefile /home/uhlemannlab/heekuk_path/Go_autocycler.smk \
#   --config Fastq_DIRS="$dir" output="$output_dir" database="$DB" threads=4 \
#   --cores 8 --rerun-incomplete --keep-going
#############################################

import os, glob


# ----------------
# CONFIG
# ----------------
FASTQ_DIRS = config.get("Fastq_DIRS", "2_qc_fastqs")
OUTDIR     = config.get("output", "3_assembled_all")

THREADS_PER_JOB = int(config.get("threads", 4))  # autocycler helper threads
JOBS            = int(config.get("jobs", 2))     # GNU parallel jobs inside autocycler
READ_TYPE       = config.get("read_type", "ont_r10")
MAX_TIME        = config.get("max_time", "12h")

QC_THREADS      = int(config.get("qc_threads", 8))
QC_MIN_Q        = int(config.get("qc_min_q", 10))
QC_MIN_LEN      = int(config.get("qc_min_len", 1000))

MEDAKA_THREADS  = int(config.get("medaka_threads", 4))

QUAST_THREADS   = int(config.get("quast_threads", 4))
QUAST_MIN_CONTIG= int(config.get("quast_min_contig", 500))

CHECKM2_THREADS = int(config.get("checkm2_threads", 8))

COV_THREADS     = int(config.get("cov_threads", 8))

BAKTA_THREADS   = int(config.get("bakta_threads", 8))

TMPDIR          = config.get("tmpdir", "")

# Single DB root (contains both bakta + checkm2 DBs)
DB_ROOT = config.get("database", "")
# porechop 조건 가동
DO_PORECHOP = int(config.get("do_porechop", 1))


# ----------------
# PATHS
# ----------------
QC_DIR         = os.path.join(OUTDIR, "1_QC")
QUAST_DIR      = os.path.join(OUTDIR, "2_quast")
AUTOCYCLER_DIR = os.path.join(OUTDIR, "3_autocycler")
MEDAKA_DIR     = os.path.join(OUTDIR, "4_medaka")
CHECKM2_DIR    = os.path.join(OUTDIR, "5_checkm2")
COV_DIR        = os.path.join(OUTDIR, "6_coverage")
BAKTA_DIR      = os.path.join(OUTDIR, "7_bakta")
SUMMARY_XLSX   = os.path.join(OUTDIR, "checkm2_coverage_summary.xlsx")


# ----------------
# DB PATHS (auto)
# ----------------
BAKTA_DB   = os.path.join(DB_ROOT, "bakta_DB") if DB_ROOT else ""
CHECKM2_DB = os.path.join(DB_ROOT, "CheckM2_database", "uniref100.KO.1.dmnd") if DB_ROOT else ""

# Safety checks (fail early)
if DB_ROOT:
    if not os.path.isdir(DB_ROOT):
        raise ValueError(f"[Go_autocycler.smk] database path does not exist: {DB_ROOT}")

    if not os.path.isdir(BAKTA_DB):
        raise ValueError(f"[Go_autocycler.smk] Bakta DB not found: {BAKTA_DB}")

    if not os.path.isfile(CHECKM2_DB):
        raise ValueError(f"[Go_autocycler.smk] CheckM2 DB not found: {CHECKM2_DB}")


# ----------------
# FASTQ discovery (raw input)
# ----------------
def discover_fastqs(fastq_dirs):
    fastq_dirs = str(fastq_dirs)
    pats = [
        os.path.join(fastq_dirs, "*.fastq.gz"),
        os.path.join(fastq_dirs, "*.fq.gz")
    ]
    files = []
    for p in pats:
        files.extend(glob.glob(p))

    def sample_name(f):
        b = os.path.basename(f)
        if b.endswith(".fastq.gz"): b = b[:-9]
        elif b.endswith(".fq.gz"): b = b[:-6]
        if b.endswith(".clean"): b = b[:-6]
        return b

    def is_clean(x):
        return x.endswith(".clean.fastq.gz") or x.endswith(".clean.fq.gz")

    m = {}
    for f in sorted(files):
        s = sample_name(f)
        if s not in m:
            m[s] = f
        else:
            if is_clean(f) and not is_clean(m[s]):
                m[s] = f
    return m


RAW_FASTQS = discover_fastqs(FASTQ_DIRS)
SAMPLES = sorted(RAW_FASTQS.keys())

if not SAMPLES:
    raise ValueError(f"[Go_autocycler.smk] No FASTQs found in {FASTQ_DIRS}")


# ----------------
# RULE ALL
# ----------------
rule all:
    input:
        # final fasta
        expand(os.path.join(MEDAKA_DIR, "{sample}_final_assembly.fasta"), sample=SAMPLES),

        # quast report
        expand(os.path.join(QUAST_DIR, "{sample}", "report.txt"), sample=SAMPLES),

        # checkm2 fixed report + done marker
        expand(os.path.join(CHECKM2_DIR, "{sample}", "quality_report.tsv"), sample=SAMPLES),
        expand(os.path.join(CHECKM2_DIR, "{sample}", "DONE.txt"), sample=SAMPLES),

        # coverage marker
        expand(os.path.join(COV_DIR, "{sample}", "DONE.txt"), sample=SAMPLES),

        # bakta marker
        expand(os.path.join(BAKTA_DIR, "{sample}", "DONE.txt"), sample=SAMPLES),

        SUMMARY_XLSX

# ----------------
# 1) QC: porechop_abi + NanoFilt
# ----------------
rule qc_clean_reads:
    input:
        fq=lambda wc: RAW_FASTQS[wc.sample]
    output:
        clean=os.path.join(QC_DIR, "{sample}.clean.fastq.gz")
    threads: QC_THREADS
    resources:
        mem_mb=8000
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{QC_DIR}"

        in_fq=$(realpath {input.fq})

        if [ -n "{TMPDIR}" ]; then
            mkdir -p "{TMPDIR}"
            export TMPDIR="{TMPDIR}"
        fi

        if [ "{DO_PORECHOP}" -eq 1 ]; then
            porechop_abi -i "$in_fq" -t {threads} \
              | NanoFilt -q {QC_MIN_Q} -l {QC_MIN_LEN} \
              | gzip > {output.clean}
        else
            zcat "$in_fq" \
              | NanoFilt -q {QC_MIN_Q} -l {QC_MIN_LEN} \
              | gzip > {output.clean}
        fi
        """


# ----------------
# 3) Autocycler assembly
# ----------------
rule autocycler_assembly:
    input:
        fq=os.path.join(QC_DIR, "{sample}.clean.fastq.gz")
    output:
        asm=os.path.join(AUTOCYCLER_DIR, "{sample}", "autocycler_out", "consensus_assembly.fasta"),
        gfa=os.path.join(AUTOCYCLER_DIR, "{sample}", "autocycler_out", "consensus_assembly.gfa")
    threads: THREADS_PER_JOB
    resources:
        mem_mb=32000
    shell:
        r"""
        set -euo pipefail
        fq_abs=$(realpath "{input.fq}")

        sample="{wildcards.sample}"
        outdir="{AUTOCYCLER_DIR}/{wildcards.sample}"
        mkdir -p "$outdir"
        cd "$outdir"

        
        fq_local=$(basename "$fq_abs")
        ln -sf "$fq_abs" "$fq_local"

        if [ -n "{TMPDIR}" ]; then
            mkdir -p "{TMPDIR}"
            export TMPDIR="{TMPDIR}"
        fi

        reads="$fq_local"
        threads="{threads}"
        jobs="{JOBS}"
        read_type="{READ_TYPE}"
        max_time="{MAX_TIME}"

        genome_size=$(autocycler helper genome_size --reads "$reads" --threads "$threads")

        autocycler subsample \
          --reads "$reads" \
          --out_dir subsampled_reads \
          --genome_size "$genome_size" \
          2>> autocycler.stderr

        mkdir -p assemblies
        rm -f assemblies/jobs.txt

        for assembler in raven miniasm flye plassembler; do
            for i in 01 02 03 04; do
                echo "autocycler helper $assembler \
                  --reads subsampled_reads/sample_$i.fastq \
                  --out_prefix assemblies/${{assembler}}_$i \
                  --threads $threads \
                  --genome_size $genome_size \
                  --read_type $read_type \
                  --min_depth_rel 0.1" >> assemblies/jobs.txt
            done
        done

        set +e
        nice -n 19 parallel \
          --jobs "$jobs" \
          --joblog assemblies/joblog.tsv \
          --results assemblies/logs \
          --timeout "$max_time" \
          < assemblies/jobs.txt
        set -e

        shopt -s nullglob
        for f in assemblies/plassembler*.fasta; do
            sed -i 's/circular=True/circular=True Autocycler_cluster_weight=3/' "$f"
        done
        for f in assemblies/flye*.fasta; do
            sed -i 's/^>.*$/& Autocycler_consensus_weight=2/' "$f"
        done
        shopt -u nullglob

        rm -f subsampled_reads/*.fastq

        autocycler compress -i assemblies -a autocycler_out 2>> autocycler.stderr
        autocycler cluster -a autocycler_out 2>> autocycler.stderr

        for c in autocycler_out/clustering/qc_pass/cluster_*; do
            autocycler trim -c "$c" 2>> autocycler.stderr
            autocycler resolve -c "$c" 2>> autocycler.stderr
        done

        autocycler combine \
          -a autocycler_out \
          -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa \
          2>> autocycler.stderr

        test -f autocycler_out/consensus_assembly.fasta
        test -f autocycler_out/consensus_assembly.gfa
        """



# ----------------
# 4) Medaka polishing
# ----------------
rule medaka_polish:
    input:
        fq=os.path.join(QC_DIR, "{sample}.clean.fastq.gz"),
        draft=os.path.join(AUTOCYCLER_DIR, "{sample}", "autocycler_out", "consensus_assembly.fasta")
    output:
        final=os.path.join(MEDAKA_DIR, "{sample}_final_assembly.fasta")
    threads: MEDAKA_THREADS
    resources:
        mem_mb=32000
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{MEDAKA_DIR}"

        fq=$(realpath {input.fq})
        draft=$(realpath {input.draft})
        out="{MEDAKA_DIR}/{wildcards.sample}_medaka"

        if [ -n "{TMPDIR}" ]; then
            mkdir -p "{TMPDIR}"
            export TMPDIR="{TMPDIR}"
        fi

        rm -rf "$out"
        mkdir -p "$out"

        medaka_consensus \
          -i "$fq" \
          -d "$draft" \
          -o "$out" \
          -t {threads} \
          --bacteria

        cp "$out/consensus.fasta" {output.final}
        """



# ----------------
# 2) QUAST (after final fasta)
# ----------------
rule quast_qc:
    input:
        asm=os.path.join(MEDAKA_DIR, "{sample}_final_assembly.fasta")
    output:
        report=os.path.join(QUAST_DIR, "{sample}", "report.txt")
    threads: QUAST_THREADS
    resources:
        mem_mb=16000
    shell:
        r"""
        set -euo pipefail
        outdir="{QUAST_DIR}/{wildcards.sample}"
        rm -rf "$outdir"
        mkdir -p "$outdir"

        micromamba run -n quast_env quast.py {input.asm} \
          -o "$outdir" \
          --threads {threads} \
          --min-contig {QUAST_MIN_CONTIG}

        test -f {output.report}
        """



# ----------------
# 5) CheckM2 (STABLE)
#   - Always produce:
#       5_checkm2/{sample}/quality_report.tsv
#       5_checkm2/{sample}/DONE.txt
#       5_checkm2/{sample}/raw/...
#
#   - Uses DB_ROOT/CheckM2_database automatically
#   - Runs via conda run -n checkm2_env (python conflict safe)
# ----------------
rule checkm2:
    input:
        asm=os.path.join(MEDAKA_DIR, "{sample}_final_assembly.fasta")
    output:
        report=os.path.join(CHECKM2_DIR, "{sample}", "quality_report.tsv"),
        done=os.path.join(CHECKM2_DIR, "{sample}", "DONE.txt")
    threads: CHECKM2_THREADS
    resources:
        mem_mb=32000
    shell:
        r"""
        set -euo pipefail

        sample="{wildcards.sample}"
        outdir="{CHECKM2_DIR}/{wildcards.sample}"
        tmpdir="$outdir/_tmp_checkm2"
        rawdir="$outdir/raw"
        logfile="$outdir/checkm2.log"

        mkdir -p "$outdir"

        # Already done? skip
        if [ -f "{output.done}" ] && [ -f "{output.report}" ]; then
            echo "[checkm2] DONE already exists for $sample, skipping."
            exit 0
        fi

        rm -rf "$tmpdir" "$rawdir"
        mkdir -p "$tmpdir"

        asm=$(realpath "{input.asm}")

        # Run checkm2
        set +e
        micromamba run -n checkm2_110 checkm2 predict \
            --input "$asm" \
            --output-directory "$tmpdir" \
            --threads {threads} \
            --database_path "{CHECKM2_DB}" \
            --force \
            > "$logfile" 2>&1
        rc=$?
        set -e

        if [ $rc -ne 0 ]; then
            echo "ERROR: CheckM2 failed for $sample (exit code=$rc)" 1>&2
            echo "----- checkm2.log -----" 1>&2
            tail -n 80 "$logfile" 1>&2 || true
            exit 1
        fi

        # Find report robustly
        report_found=$(find "$tmpdir" -maxdepth 4 -type f -name "quality_report.tsv" | head -n 1 || true)

        if [ -z "$report_found" ]; then
            report_found=$(find "$tmpdir" -maxdepth 4 -type f -name "*quality_report*.tsv" | head -n 1 || true)
        fi

        if [ -z "$report_found" ]; then
            echo "ERROR: CheckM2 finished but no quality_report.tsv found for $sample" 1>&2
            echo "Files under tmpdir:" 1>&2
            find "$tmpdir" -maxdepth 4 -type f 1>&2 || true
            echo "----- checkm2.log -----" 1>&2
            tail -n 80 "$logfile" 1>&2 || true
            exit 1
        fi

        cp "$report_found" "{output.report}"
        mv "$tmpdir" "$rawdir"

        echo "DONE" > "{output.done}"
        """



# ----------------
# 6) Coverage (minimap2 + samtools + contig mean depth)
# ----------------
rule coverage:
    input:
        fq=os.path.join(QC_DIR, "{sample}.clean.fastq.gz"),
        asm=os.path.join(MEDAKA_DIR, "{sample}_final_assembly.fasta")
    output:
        bam=os.path.join(COV_DIR, "{sample}", "{sample}.sorted.bam"),
        bai=os.path.join(COV_DIR, "{sample}", "{sample}.sorted.bam.bai"),
        tsv=os.path.join(COV_DIR, "{sample}", "{sample}_contig_mean_depth.tsv"),
        done=os.path.join(COV_DIR, "{sample}", "DONE.txt")
    threads: COV_THREADS
    resources:
        mem_mb=32000
    shell:
        r"""
        set -euo pipefail
        outdir="{COV_DIR}/{wildcards.sample}"
        rm -rf "$outdir"
        mkdir -p "$outdir"

        fq=$(realpath {input.fq})
        asm=$(realpath {input.asm})

        minimap2 -ax map-ont -t {threads} "$asm" "$fq" \
          | samtools sort -@ {threads} -o "{output.bam}" -

        samtools index "{output.bam}"

        samtools depth -a "{output.bam}" \
          | awk '{{sum[$1]+=$3; cnt[$1]++}} END {{for (c in sum) printf "%s\t%.3f\n", c, sum[c]/cnt[c]}}' \
          | sort -k2,2nr > "{output.tsv}"

        echo "DONE" > "{output.done}"
        """



# ----------------
# 7) Bakta
# ----------------
rule bakta_annotate:
    input:
        fasta=os.path.join(MEDAKA_DIR, "{sample}_final_assembly.fasta")
    output:
        done=os.path.join(BAKTA_DIR, "{sample}", "DONE.txt")
    threads: BAKTA_THREADS
    resources:
        mem_mb=32000
    shell:
        r"""
        set -euo pipefail
        sample="{wildcards.sample}"
        outdir="{BAKTA_DIR}/{wildcards.sample}"

        mkdir -p "{BAKTA_DIR}"

        fasta=$(realpath "{input.fasta}")

        bakta \
          --db "{BAKTA_DB}" \
          --output "$outdir" \
          --prefix "$sample" \
          --threads {threads} \
          --force \
          "$fasta"

        echo "DONE" > "{output.done}"
        """

        
# ----------------
# 9) Summary Excel (CheckM2 + Coverage)
# ----------------
SUMMARY_XLSX = os.path.join(OUTDIR, "checkm2_coverage_summary.xlsx")

rule make_checkm2_coverage_summary:
    input:
        checkm2=expand(os.path.join(CHECKM2_DIR, "{sample}", "quality_report.tsv"), sample=SAMPLES),
        cov=expand(os.path.join(COV_DIR, "{sample}", "{sample}_contig_mean_depth.tsv"), sample=SAMPLES)
    output:
        xlsx=SUMMARY_XLSX
    threads: 1
    resources:
        mem_mb=4000
    shell:
        r"""
        set -euo pipefail

        python - << 'PY'
import os
import glob
import pandas as pd

checkm2_dir = "{CHECKM2_DIR}"
cov_dir = "{COV_DIR}"
out_xlsx = "{output.xlsx}"

# -------------------------
# 1) CheckM2 merge
# -------------------------
checkm2_files = sorted(glob.glob(os.path.join(checkm2_dir, "*", "quality_report.tsv")))

checkm2_all = []
for f in checkm2_files:
    sample = os.path.basename(os.path.dirname(f))
    df = pd.read_csv(f, sep="\t")
    df.insert(0, "sample", sample)
    checkm2_all.append(df)

checkm2_all = pd.concat(checkm2_all, ignore_index=True) if checkm2_all else pd.DataFrame()

# -------------------------
# 2) Coverage summary
# -------------------------
cov_files = sorted(glob.glob(os.path.join(cov_dir, "*", "*_contig_mean_depth.tsv")))

cov_rows = []
for f in cov_files:
    sample = os.path.basename(os.path.dirname(f))
    df = pd.read_csv(f, sep="\t", header=None, names=["contig", "mean_depth"])

    if df.shape[0] == 0:
        cov_rows.append({{"sample": sample}})
        continue

    df_sorted = df.sort_values("mean_depth", ascending=False)

    cov_rows.append({{
        "sample": sample,
        "n_contigs": int(df.shape[0]),
        "mean_depth_overall": float(df["mean_depth"].mean()),
        "median_depth_overall": float(df["mean_depth"].median()),
        "max_depth": float(df["mean_depth"].max()),
        "min_depth": float(df["mean_depth"].min()),
        "top_contig": str(df_sorted.iloc[0]["contig"]),
        "top_contig_depth": float(df_sorted.iloc[0]["mean_depth"]),
    }})

coverage_summary = pd.DataFrame(cov_rows)

# -------------------------
# 3) Write Excel
# -------------------------
outdir = os.path.dirname(out_xlsx)
if outdir:
    os.makedirs(outdir, exist_ok=True)

with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
    checkm2_all.to_excel(writer, sheet_name="CheckM2", index=False)
    coverage_summary.to_excel(writer, sheet_name="Coverage", index=False)

print("[OK] wrote:", out_xlsx)
PY
        """
