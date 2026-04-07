#############################################
# Go_longWGS_V6_3_docker.smk
# 20260401
# Heekuk Park
#
# Pipeline:
#   QC → autocycler → medaka → quast → checkm2 → coverage → kraken2 → bakta → bandage
#
# OUTDIR/
# └── 1_QC/             (porechop_abi + NanoFilt output: *.clean.fastq.gz)
# └── 2_quast/          (QUAST reports; from medaka fasta)
# └── 3_autocycler/     ({sample}/autocycler_out/consensus_assembly.{fasta,gfa})
# └── 4_medaka/         ({sample}_final_assembly.fasta)
# └── 5_checkm2/        ({sample}/quality_report.tsv + DONE.txt + raw/)
# └── 6_coverage/       ({sample}.sorted.bam + {sample}_contig_mean_depth.tsv + DONE.txt)
# └── 7a_kraken2/       ({sample}.kreport + kraken2_summary.csv + db_info)
# └── 7_bakta/          ({sample}/DONE.txt + outputs)
#
# DB structure (single root):
#   DB_ROOT/
#     ├── bakta_DB/
#     └── CheckM2_database/
#     └── kraken2DB/ or kraken2_db/ or kraken_db/
#         └── <kraken_db_files>
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

import os, glob, json, shlex


# ----------------
# CONFIG
# ----------------
FASTQ_DIRS = config.get("Fastq_DIRS", "2_qc_fastqs")
OUTDIR     = config.get("output", "3_assembled_all")
PIPELINE_MODE = str(config.get("pipeline_mode", "strict")).strip().lower()

THREADS_PER_JOB = int(config.get("threads", 2))  # autocycler helper threads
JOBS            = int(config.get("jobs", 1))     # GNU parallel jobs inside autocycler
READ_TYPE       = config.get("read_type", "ont_r10")
MAX_TIME        = config.get("max_time", "12h")

QC_THREADS      = int(config.get("qc_threads", 2))
QC_MIN_Q        = int(config.get("qc_min_q", 10))
QC_MIN_LEN      = int(config.get("qc_min_len", 1000))

MEDAKA_THREADS  = int(config.get("medaka_threads", 2))
MEDAKA_MODEL    = str(config.get("medaka_model", "")).strip()
MEDAKA_FALLBACK_MODEL = str(
    config.get("medaka_fallback_model", MEDAKA_MODEL if MEDAKA_MODEL else "r1041_e82_400bps_hac_g632")
).strip()
DNAAPLER_THREADS = int(config.get("dnaapler_threads", MEDAKA_THREADS))

QUAST_THREADS   = int(config.get("quast_threads", 2))
QUAST_MIN_CONTIG= int(config.get("quast_min_contig", 500))

CHECKM2_THREADS = int(config.get("checkm2_threads", 2))

COV_THREADS     = int(config.get("cov_threads", 2))

KRAKEN_THREADS  = int(config.get("kraken_threads", 2))

BAKTA_THREADS   = int(config.get("bakta_threads", 2))

# Stage memory limits (override via --config)
QC_MEM_MB         = int(config.get("qc_mem_mb", 8000))
AUTOCYCLER_MEM_MB = int(config.get("autocycler_mem_mb", 32000))
MEDAKA_MEM_MB     = int(config.get("medaka_mem_mb", 32000))
DNAAPLER_MEM_MB   = int(config.get("dnaapler_mem_mb", 8000))
QUAST_MEM_MB      = int(config.get("quast_mem_mb", 16000))
CHECKM2_MEM_MB    = int(config.get("checkm2_mem_mb", 32000))
COV_MEM_MB        = int(config.get("cov_mem_mb", 32000))
KRAKEN_MEM_MB     = int(config.get("kraken_mem_mb", 32000))
BAKTA_MEM_MB      = int(config.get("bakta_mem_mb", 32000))
SUMMARY_MEM_MB    = int(config.get("summary_mem_mb", 4000))

TMPDIR          = config.get("tmpdir", "")

# Single DB root (contains bakta + checkm2 + optional kraken2 DB)
DB_ROOT = config.get("database", "")
# porechop 조건 가동
DO_PORECHOP = int(config.get("do_porechop", 1))


def _is_kraken_db_dir(path):
    if not path or not os.path.isdir(path):
        return False
    # Kraken2 DB layouts vary a bit depending on build/distribution style.
    # Accept either the classic trio (opts/hash/taxo) or lighter layouts that
    # still contain the core hash file plus taxonomy/index sidecars.
    markers = {name: os.path.exists(os.path.join(path, name)) for name in [
        "hash.k2d",
        "opts.k2d",
        "taxo.k2d",
        "nodes.dmp",
        "seqid2taxid.map",
        "inspect.txt",
    ]}

    if markers["hash.k2d"] and (markers["opts.k2d"] or markers["taxo.k2d"]):
        return True

    if markers["hash.k2d"] and (markers["nodes.dmp"] or markers["seqid2taxid.map"] or markers["inspect.txt"]):
        return True

    return False


def _find_kraken_db(db_root):
    explicit = str(config.get("kraken_db", "")).strip()
    if explicit:
        return explicit

    if not db_root:
        return ""

    candidates = [
        os.path.join(db_root, "kraken2DB"),
        os.path.join(db_root, "kraken2_db"),
        os.path.join(db_root, "kraken_db"),
        os.path.join(db_root, "Kraken2DB"),
        os.path.join(db_root, "Kraken2_db"),
    ]

    for cand in candidates:
        if _is_kraken_db_dir(cand):
            return cand

    child_hits = []
    for cand in candidates:
        if os.path.isdir(cand):
            for child in sorted(glob.glob(os.path.join(cand, "*"))):
                if _is_kraken_db_dir(child):
                    child_hits.append(child)
    if len(child_hits) == 1:
        return child_hits[0]

    shallow_hits = []
    for child in sorted(glob.glob(os.path.join(db_root, "*"))):
        if _is_kraken_db_dir(child):
            shallow_hits.append(child)
    if len(shallow_hits) == 1:
        return shallow_hits[0]

    return ""


# ----------------
# PATHS
# ----------------
QC_DIR         = os.path.join(OUTDIR, "1_QC")
QUAST_DIR      = os.path.join(OUTDIR, "2_quast")
AUTOCYCLER_DIR = os.path.join(OUTDIR, "3_autocycler")
MEDAKA_DIR     = os.path.join(OUTDIR, "4_medaka")
CHECKM2_DIR    = os.path.join(OUTDIR, "5_checkm2")
COV_DIR        = os.path.join(OUTDIR, "6_coverage")
KRAKEN_DIR     = os.path.join(OUTDIR, "7a_kraken2")
BAKTA_DIR      = os.path.join(OUTDIR, "7_bakta")
SUMMARY_XLSX   = os.path.join(OUTDIR, "checkm2_coverage_summary.xlsx")
AUTOCYCLER_TABLE_TSV = os.path.join(OUTDIR, "autocycler_metrics.tsv")
FAIL_LOG       = os.path.join(OUTDIR, "0_failed_samples.tsv")


# ----------------
# DB PATHS (auto)
# ----------------
BAKTA_DB   = os.path.join(DB_ROOT, "bakta_DB") if DB_ROOT else ""
CHECKM2_DB = os.path.join(DB_ROOT, "CheckM2_database", "uniref100.KO.1.dmnd") if DB_ROOT else ""
KRAKEN_DB  = _find_kraken_db(DB_ROOT)
KRAKEN_DB_NAME = os.path.basename(KRAKEN_DB.rstrip("/")) if KRAKEN_DB else ""

# Safety checks (fail early)
if DB_ROOT:
    if not os.path.isdir(DB_ROOT):
        raise ValueError(f"[Go_autocycler.smk] database path does not exist: {DB_ROOT}")

    if not os.path.isdir(BAKTA_DB):
        raise ValueError(f"[Go_autocycler.smk] Bakta DB not found: {BAKTA_DB}")

    if not os.path.isfile(CHECKM2_DB):
        raise ValueError(f"[Go_autocycler.smk] CheckM2 DB not found: {CHECKM2_DB}")

    if not KRAKEN_DB:
        raise ValueError(
            "[Go_autocycler.smk] Kraken2 DB not found under DB root. "
            "Expected e.g. DB_ROOT/kraken2DB/<db> or set --config kraken_db=/db/<path>"
        )


# ----------------
# FASTQ discovery (raw input)
# ----------------
def discover_fastqs(fastq_dirs):
    if isinstance(fastq_dirs, (list, tuple)):
        dirs = [str(x).strip() for x in fastq_dirs if str(x).strip()]
    else:
        dirs = [x.strip() for x in str(fastq_dirs).split(",") if x.strip()]

    pats = []
    for d in dirs:
        pats.extend([
            os.path.join(d, "*.fastq.gz"),
            os.path.join(d, "*.fq.gz")
        ])
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
STRICT_ALL_INPUTS = [
    # final fasta
    *expand(os.path.join(MEDAKA_DIR, "{sample}_final_assembly.fasta"), sample=SAMPLES),
    # quast report
    *expand(os.path.join(QUAST_DIR, "{sample}", "report.txt"), sample=SAMPLES),
    # checkm2 fixed report + done marker
    *expand(os.path.join(CHECKM2_DIR, "{sample}", "quality_report.tsv"), sample=SAMPLES),
    *expand(os.path.join(CHECKM2_DIR, "{sample}", "DONE.txt"), sample=SAMPLES),
    # coverage marker
    *expand(os.path.join(COV_DIR, "{sample}", "DONE.txt"), sample=SAMPLES),
    # kraken2 outputs
    *expand(os.path.join(KRAKEN_DIR, "{sample}.done.txt"), sample=SAMPLES),
    os.path.join(KRAKEN_DIR, "kraken2_summary.csv"),
    os.path.join(KRAKEN_DIR, "kraken2_db_info.txt"),
    # bakta marker
    *expand(os.path.join(BAKTA_DIR, "{sample}", "DONE.txt"), sample=SAMPLES),
    AUTOCYCLER_TABLE_TSV,
    SUMMARY_XLSX
]

PERMISSIVE_ALL_INPUTS = [AUTOCYCLER_TABLE_TSV, SUMMARY_XLSX]

ALL_INPUTS = PERMISSIVE_ALL_INPUTS if PIPELINE_MODE == "permissive" else STRICT_ALL_INPUTS

rule all:
    input:
        ALL_INPUTS

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
        mem_mb=QC_MEM_MB
    shell:
        r"""
        set -euo pipefail
        fail_log="{FAIL_LOG}"
        mkdir -p "$(dirname "$fail_log")"
        [ -f "$fail_log" ] || echo -e "sample\tstage\treason" > "$fail_log"
        trap 'echo -e "{wildcards.sample}\tqc\tcommand_failed" >> "$fail_log"' ERR
        mkdir -p "{QC_DIR}"

        in_fq=$(realpath {input.fq})
        gzip -t "$in_fq"

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
        mem_mb=AUTOCYCLER_MEM_MB
    shell:
        r"""
        set -euo pipefail
        fail_log=$(realpath -m "{FAIL_LOG}")
        mkdir -p "$(dirname "$fail_log")"
        [ -f "$fail_log" ] || echo -e "sample\tstage\treason" > "$fail_log"
        echo "[autocycler] START sample={wildcards.sample}"
        fq_abs=$(realpath "{input.fq}")

        outdir_abs=$(realpath -m "{AUTOCYCLER_DIR}/{wildcards.sample}")
        mkdir -p "$outdir_abs"
        cd "$outdir_abs"

        fail_and_exit() {{
            local reason="$1"
            local size_mb="NA"
            local last_err="NA"
            trap - ERR
            set +e
            if [ -n "${{fq_abs:-}}" ] && [ -f "$fq_abs" ]; then
                size_mb=$(du -m "$fq_abs" 2>/dev/null | awk '{{print $1}}')
            fi
            if [ -f autocycler.stderr ]; then
                last_err=$(tail -n 1 autocycler.stderr | tr '\t' ' ' | tr '\r\n' ' ')
                [ -n "$last_err" ] || last_err="NA"
            fi
            echo -e "{wildcards.sample}\tautocycler\t$reason;size_mb=${{size_mb}};last_err=${{last_err}}" >> "$fail_log" || true
            echo "[autocycler][ERROR] sample={wildcards.sample} failed ($reason)." 1>&2
            tail -n 120 autocycler.stderr 1>&2 || true
            exit 1
        }}
        trap 'fail_and_exit "command_failed"' ERR

        
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

        genome_size=$(autocycler helper genome_size --reads "$reads" --threads "$threads" 2>> autocycler.stderr)
        echo "[autocycler] genome_size=$genome_size"

        echo "[autocycler] subsample ..."
        autocycler subsample \
          --reads "$reads" \
          --out_dir subsampled_reads \
          --genome_size "$genome_size" \
          2>> autocycler.stderr || fail_and_exit "subsample_failed"

        echo "[autocycler] prepare full-read plassembler input ..."
        full_reads_fastq="full_reads_for_plassembler.fastq"
        zcat "$reads" > "$full_reads_fastq"

        mkdir -p assemblies
        rm -f assemblies/jobs.txt

        for assembler in raven miniasm flye; do
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

        # Keep plasmid-sensitive assembly on the full QC read set so low-abundance
        # plasmid reads are not lost during the Autocycler subsampling step.
        echo "autocycler helper plassembler \
          --reads $full_reads_fastq \
          --out_prefix assemblies/plassembler_full \
          --threads $threads \
          --genome_size $genome_size \
          --read_type $read_type \
          --min_depth_rel 0.1" >> assemblies/jobs.txt

        set +e
        nice -n 19 parallel \
          --jobs "$jobs" \
          --joblog assemblies/joblog.tsv \
          --results assemblies/logs \
          --timeout "$max_time" \
          < assemblies/jobs.txt
        prc=$?
        set -e
        if [ "$prc" -ne 0 ]; then
            echo "[autocycler][ERROR] parallel returned non-zero: $prc" 1>&2
            if [ -f assemblies/joblog.tsv ]; then
                fail_n=$(awk -F'\t' 'NR>1 && $7 != 0 {{n++}} END{{print n+0}}' assemblies/joblog.tsv)
                echo "[autocycler][ERROR] failed helper jobs: $fail_n" 1>&2
                awk -F'\t' 'NR==1 || (NR>1 && $7 != 0)' assemblies/joblog.tsv | head -n 20 1>&2 || true
            fi
            tail -n 60 autocycler.stderr 1>&2 || true
            fail_and_exit "parallel_nonzero_$prc"
        fi

        shopt -s nullglob
        for f in assemblies/plassembler*.fasta; do
            sed -i 's/circular=True/circular=True Autocycler_cluster_weight=3/' "$f"
        done
        for f in assemblies/flye*.fasta; do
            sed -i 's/^>.*$/& Autocycler_consensus_weight=2/' "$f"
        done
        shopt -u nullglob

        rm -f subsampled_reads/*.fastq "$full_reads_fastq"

        echo "[autocycler] compress ..."
        autocycler compress -i assemblies -a autocycler_out 2>> autocycler.stderr || fail_and_exit "compress_failed"
        echo "[autocycler] cluster ..."
        autocycler cluster -a autocycler_out 2>> autocycler.stderr || fail_and_exit "cluster_failed"

        shopt -s nullglob
        clusters=(autocycler_out/clustering/qc_pass/cluster_*)
        shopt -u nullglob
        if [ "${{#clusters[@]}}" -eq 0 ]; then
            echo "[autocycler][ERROR] No qc_pass clusters found. Possible fragmented/contaminated assemblies." 1>&2
            tail -n 120 autocycler.stderr 1>&2 || true
            fail_and_exit "no_qc_pass_clusters"
        fi

        echo "[autocycler] trim/resolve clusters: ${{#clusters[@]}}"
        for c in "${{clusters[@]}}"; do
            autocycler trim -c "$c" 2>> autocycler.stderr || fail_and_exit "trim_failed"
            autocycler resolve -c "$c" 2>> autocycler.stderr || fail_and_exit "resolve_failed"
        done

        echo "[autocycler] combine ..."
        autocycler combine \
          -a autocycler_out \
          -i "${{clusters[@]/%//5_final.gfa}}" \
          2>> autocycler.stderr || fail_and_exit "combine_failed"

        test -f autocycler_out/consensus_assembly.fasta || fail_and_exit "missing_consensus_fasta"
        test -f autocycler_out/consensus_assembly.gfa || fail_and_exit "missing_consensus_gfa"
        echo "[autocycler] DONE sample={wildcards.sample}"
        """


rule autocycler_reorient_with_dnaapler:
    input:
        fasta=os.path.join(AUTOCYCLER_DIR, "{sample}", "autocycler_out", "consensus_assembly.fasta"),
        gfa=os.path.join(AUTOCYCLER_DIR, "{sample}", "autocycler_out", "consensus_assembly.gfa")
    output:
        fasta=os.path.join(AUTOCYCLER_DIR, "{sample}", "autocycler_out", "consensus_assembly.reoriented.fasta"),
        done=os.path.join(AUTOCYCLER_DIR, "{sample}", "autocycler_out", "dnaapler.done.txt")
    threads: DNAAPLER_THREADS
    resources:
        mem_mb=DNAAPLER_MEM_MB
    shell:
        r"""
        set -euo pipefail
        fail_log=$(realpath -m "{FAIL_LOG}")
        mkdir -p "$(dirname "$fail_log")"
        [ -f "$fail_log" ] || echo -e "sample\tstage\treason" > "$fail_log"
        trap 'echo -e "{wildcards.sample}\tdnaapler\tcommand_failed" >> "$fail_log"' ERR

        in_fasta=$(realpath "{input.fasta}")
        outdir=$(dirname "{output.fasta}")
        tmpdir="$outdir/dnaapler_out"
        mkdir -p "$outdir"
        rm -rf "$tmpdir"
        rm -f "{output.fasta}" "{output.done}"

        dnaapler all \
          -i "$in_fasta" \
          -o "$tmpdir" \
          -p "{wildcards.sample}" \
          -t {threads}

        reoriented_fasta=$(find "$tmpdir" -type f -name "*.fasta" ! -path "*/log/*" | head -n 1 || true)
        if [ -z "$reoriented_fasta" ]; then
            echo -e "{wildcards.sample}\tdnaapler\treoriented_fasta_missing" >> "$fail_log"
            echo "[dnaapler][ERROR] Reoriented FASTA not found for {wildcards.sample}" 1>&2
            find "$tmpdir" -maxdepth 3 -type f 1>&2 || true
            exit 1
        fi

        cp "$reoriented_fasta" "{output.fasta}"
        echo "DONE" > "{output.done}"
        """


rule autocycler_make_metrics_table:
    input:
        expand(os.path.join(AUTOCYCLER_DIR, "{sample}", "autocycler_out", "consensus_assembly.gfa"), sample=SAMPLES)
    params:
        sample_dirs=" ".join(shlex.quote(os.path.join(AUTOCYCLER_DIR, sample)) for sample in SAMPLES)
    output:
        tsv=AUTOCYCLER_TABLE_TSV
    threads: 1
    resources:
        mem_mb=SUMMARY_MEM_MB
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.tsv}")"
        raw_tsv="{output.tsv}.raw"
        : > "$raw_tsv"

        for sample_dir in {params.sample_dirs}; do
            sample_name=$(basename "$sample_dir")
            autocycler table -a "$sample_dir" -n "$sample_name" >> "$raw_tsv"
        done

        python - << 'PY'
from pathlib import Path
import csv

raw_path = Path(r"{output.tsv}.raw")
out_path = Path(r"{output.tsv}")

rows = []
with raw_path.open() as handle:
    for line in handle:
        line = line.rstrip("\n")
        if not line.strip():
            continue
        rows.append(line.split("\t"))

if not rows:
    out_path.write_text("")
    raise SystemExit(0)

def looks_like_header(row):
    if len(row) < 2:
        return False
    first = row[0].strip().lower()
    second = row[1].strip().lower()
    header_tokens = {"name", "sample", "read_count", "reads", "input_read_count"}
    return first in header_tokens or second in header_tokens

if looks_like_header(rows[0]):
    header = rows[0]
    data_rows = rows[1:]
else:
    max_cols = max(len(r) for r in rows)
    header = ["sample"] + [f"metric_{i:02d}" for i in range(1, max_cols)]
    data_rows = rows

with out_path.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(header)
    for row in data_rows:
        padded = row + [""] * (len(header) - len(row))
        if padded == header:
            continue
        writer.writerow(padded[:len(header)])
PY

        rm -f "$raw_tsv"
        """



# ----------------
# 4) Medaka polishing
# ----------------
rule medaka_polish:
    input:
        fq=os.path.join(QC_DIR, "{sample}.clean.fastq.gz"),
        draft=os.path.join(AUTOCYCLER_DIR, "{sample}", "autocycler_out", "consensus_assembly.reoriented.fasta")
    output:
        final=os.path.join(MEDAKA_DIR, "{sample}_final_assembly.fasta")
    threads: MEDAKA_THREADS
    resources:
        mem_mb=MEDAKA_MEM_MB
    shell:
        r"""
        set -euo pipefail
        fail_log=$(realpath -m "{FAIL_LOG}")
        mkdir -p "$(dirname "$fail_log")"
        [ -f "$fail_log" ] || echo -e "sample\tstage\treason" > "$fail_log"
        trap 'echo -e "{wildcards.sample}\tmedaka\tcommand_failed" >> "$fail_log"' ERR
        mkdir -p "{MEDAKA_DIR}"

        fq=$(realpath {input.fq})
        draft=$(realpath {input.draft})
        out="{MEDAKA_DIR}/{wildcards.sample}_medaka"
        medaka_log="$out/medaka.log"
        fallback_model="{MEDAKA_FALLBACK_MODEL}"
        draft_local=""

        if [ -n "{TMPDIR}" ]; then
            mkdir -p "{TMPDIR}"
            export TMPDIR="{TMPDIR}"
        fi

        rm -rf "$out"
        mkdir -p "$out"
        draft_local="$out/draft.fasta"
        cp "$draft" "$draft_local"
        rm -f "$draft_local.fai" "$draft_local.map-ont.mmi"

        set +e
        medaka_consensus \
          -i "$fq" \
          -d "$draft_local" \
          -o "$out" \
          -t {threads} \
          --bacteria \
          > "$medaka_log" 2>&1
        rc=$?
        set -e

        if [ "$rc" -ne 0 ]; then
            echo "[medaka][WARN] --bacteria failed for {wildcards.sample}; retrying with -m $fallback_model" 1>&2
            rm -rf "$out"
            mkdir -p "$out"
            medaka_log="$out/medaka.log"
            draft_local="$out/draft.fasta"
            cp "$draft" "$draft_local"
            rm -f "$draft_local.fai" "$draft_local.map-ont.mmi"

            set +e
            medaka_consensus \
              -i "$fq" \
              -d "$draft_local" \
              -o "$out" \
              -t {threads} \
              -m "$fallback_model" \
              > "$medaka_log" 2>&1
            rc=$?
            set -e

            if [ "$rc" -ne 0 ]; then
                echo "[medaka][ERROR] fallback model failed for {wildcards.sample}: $fallback_model" 1>&2
                tail -n 80 "$medaka_log" 1>&2 || true
                exit 1
            fi
        fi

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
        mem_mb=QUAST_MEM_MB
    shell:
        r"""
        set -euo pipefail
        fail_log=$(realpath -m "{FAIL_LOG}")
        mkdir -p "$(dirname "$fail_log")"
        [ -f "$fail_log" ] || echo -e "sample\tstage\treason" > "$fail_log"
        trap 'echo -e "{wildcards.sample}\tquast\tcommand_failed" >> "$fail_log"' ERR
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
        mem_mb=CHECKM2_MEM_MB
    shell:
        r"""
        set -euo pipefail
        fail_log=$(realpath -m "{FAIL_LOG}")
        mkdir -p "$(dirname "$fail_log")"
        [ -f "$fail_log" ] || echo -e "sample\tstage\treason" > "$fail_log"
        trap 'echo -e "{wildcards.sample}\tcheckm2\tcommand_failed" >> "$fail_log"' ERR

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
            echo -e "{wildcards.sample}\tcheckm2\tpredict_nonzero_$rc" >> "$fail_log"
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
            echo -e "{wildcards.sample}\tcheckm2\treport_missing" >> "$fail_log"
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
        mem_mb=COV_MEM_MB
    shell:
        r"""
        set -euo pipefail
        fail_log=$(realpath -m "{FAIL_LOG}")
        mkdir -p "$(dirname "$fail_log")"
        [ -f "$fail_log" ] || echo -e "sample\tstage\treason" > "$fail_log"
        trap 'echo -e "{wildcards.sample}\tcoverage\tcommand_failed" >> "$fail_log"' ERR
        outdir="{COV_DIR}/{wildcards.sample}"
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
# 7) Kraken2 (clean long reads)
# ----------------
rule kraken2_classify:
    input:
        fq=os.path.join(QC_DIR, "{sample}.clean.fastq.gz")
    output:
        kreport=os.path.join(KRAKEN_DIR, "{sample}.kreport"),
        out=os.path.join(KRAKEN_DIR, "{sample}.kraken"),
        done=os.path.join(KRAKEN_DIR, "{sample}.done.txt")
    threads: KRAKEN_THREADS
    resources:
        mem_mb=KRAKEN_MEM_MB
    shell:
        r"""
        set -euo pipefail
        fail_log=$(realpath -m "{FAIL_LOG}")
        mkdir -p "$(dirname "$fail_log")"
        [ -f "$fail_log" ] || echo -e "sample\tstage\treason" > "$fail_log"
        trap 'echo -e "{wildcards.sample}\tkraken2\tcommand_failed" >> "$fail_log"' ERR
        mkdir -p "{KRAKEN_DIR}"

        fq=$(realpath "{input.fq}")

        kraken2 \
          --db "{KRAKEN_DB}" \
          --gzip-compressed \
          --threads {threads} \
          --report "{output.kreport}" \
          --output "{output.out}" \
          "$fq"

        echo "DONE" > "{output.done}"
        """


rule kraken2_db_info:
    output:
        info=os.path.join(KRAKEN_DIR, "kraken2_db_info.txt")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{KRAKEN_DIR}"
        printf "kraken_db_path\t%s\n" "{KRAKEN_DB}" > "{output.info}"
        printf "kraken_db_name\t%s\n" "{KRAKEN_DB_NAME}" >> "{output.info}"
        """


rule kraken2_make_summary_table:
    input:
        kreports=expand(os.path.join(KRAKEN_DIR, "{sample}.kreport"), sample=SAMPLES),
        db_info=os.path.join(KRAKEN_DIR, "kraken2_db_info.txt")
    output:
        table=os.path.join(KRAKEN_DIR, "kraken2_summary.csv")
    run:
        import csv

        os.makedirs(KRAKEN_DIR, exist_ok=True)

        def parse_kreport(path):
            top_species = []
            unclassified_pct = 0.0
            try:
                with open(path) as fh:
                    for line in fh:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 6:
                            continue
                        try:
                            pct = float(parts[0].strip())
                        except Exception:
                            continue
                        rank = parts[3].strip()
                        name = parts[5].strip()
                        if rank == "U":
                            unclassified_pct = pct
                        elif rank == "S":
                            top_species.append((pct, name))
            except Exception:
                return [], 0.0

            top_species.sort(key=lambda x: x[0], reverse=True)
            return top_species[:3], unclassified_pct

        def fmt_species(item):
            if not item:
                return "NA"
            pct, name = item
            return f"{name} ({pct:.2f}%)"

        rows = []
        for sample in SAMPLES:
            kreport = os.path.join(KRAKEN_DIR, f"{sample}.kreport")
            top3, unclassified_pct = parse_kreport(kreport)
            rows.append({
                "Sample": sample,
                "Top1": fmt_species(top3[0] if len(top3) > 0 else None),
                "Top2": fmt_species(top3[1] if len(top3) > 1 else None),
                "Top3": fmt_species(top3[2] if len(top3) > 2 else None),
                "Unclassified": f"{unclassified_pct:.2f}%"
            })

        with open(output.table, "w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=["Sample", "Top1", "Top2", "Top3", "Unclassified"])
            writer.writeheader()
            writer.writerows(rows)


# ----------------
# 8) Bakta
# ----------------
rule bakta_annotate:
    input:
        fasta=os.path.join(MEDAKA_DIR, "{sample}_final_assembly.fasta")
    output:
        done=os.path.join(BAKTA_DIR, "{sample}", "DONE.txt")
    threads: BAKTA_THREADS
    resources:
        mem_mb=BAKTA_MEM_MB
    shell:
        r"""
        set -euo pipefail
        fail_log=$(realpath -m "{FAIL_LOG}")
        mkdir -p "$(dirname "$fail_log")"
        [ -f "$fail_log" ] || echo -e "sample\tstage\treason" > "$fail_log"
        trap 'echo -e "{wildcards.sample}\tbakta\tcommand_failed" >> "$fail_log"' ERR
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
# 9) Summary Excel (CheckM2 + Coverage + Kraken2 + QUAST)
# ----------------
def existing_summary_inputs(_wc=None):
    deps = []
    for sample in SAMPLES:
        for p in [
            os.path.join(CHECKM2_DIR, sample, "DONE.txt"),
            os.path.join(COV_DIR, sample, "DONE.txt"),
            os.path.join(KRAKEN_DIR, f"{sample}.done.txt"),
            os.path.join(BAKTA_DIR, sample, "DONE.txt"),
            os.path.join(QUAST_DIR, sample, "report.txt"),
            os.path.join(MEDAKA_DIR, f"{sample}_final_assembly.fasta"),
            os.path.join(AUTOCYCLER_DIR, sample, "autocycler_out", "dnaapler.done.txt"),
        ]:
            if os.path.exists(p):
                deps.append(p)
    for p in [
        AUTOCYCLER_TABLE_TSV,
        os.path.join(KRAKEN_DIR, "kraken2_db_info.txt"),
        os.path.join(KRAKEN_DIR, "kraken2_summary.csv"),
    ]:
        if os.path.exists(p):
            deps.append(p)
    return sorted(set(deps))


_STRICT_SUMMARY_DEPS = (
    expand(os.path.join(CHECKM2_DIR, "{sample}", "DONE.txt"), sample=SAMPLES) +
    expand(os.path.join(COV_DIR, "{sample}", "DONE.txt"), sample=SAMPLES) +
    expand(os.path.join(KRAKEN_DIR, "{sample}.done.txt"), sample=SAMPLES) +
    [os.path.join(KRAKEN_DIR, "kraken2_summary.csv"),
     os.path.join(KRAKEN_DIR, "kraken2_db_info.txt")] +
    expand(os.path.join(BAKTA_DIR, "{sample}", "DONE.txt"), sample=SAMPLES) +
    expand(os.path.join(QUAST_DIR, "{sample}", "report.txt"), sample=SAMPLES) +
    expand(os.path.join(MEDAKA_DIR, "{sample}_final_assembly.fasta"), sample=SAMPLES) +
    expand(os.path.join(AUTOCYCLER_DIR, "{sample}", "autocycler_out", "dnaapler.done.txt"), sample=SAMPLES) +
    [AUTOCYCLER_TABLE_TSV]
)


rule make_checkm2_coverage_summary:
    input:
        existing_summary_inputs if PIPELINE_MODE == "permissive" else _STRICT_SUMMARY_DEPS
    params:
        samples_json=lambda wc: json.dumps(SAMPLES)
    output:
        xlsx=SUMMARY_XLSX
    threads: 1
    resources:
        mem_mb=SUMMARY_MEM_MB
    shell:
        r"""
        set -euo pipefail

        python - << 'PY'
import os
import glob
import json
import pandas as pd
from pathlib import Path

checkm2_dir = "{CHECKM2_DIR}"
cov_dir = "{COV_DIR}"
qc_dir = "{QC_DIR}"
autocycler_dir = "{AUTOCYCLER_DIR}"
medaka_dir = "{MEDAKA_DIR}"
quast_dir = "{QUAST_DIR}"
bakta_dir = "{BAKTA_DIR}"
kraken_dir = "{KRAKEN_DIR}"
autocycler_table_tsv = "{AUTOCYCLER_TABLE_TSV}"
outdir_root = "{OUTDIR}"
out_xlsx = "{output.xlsx}"
samples = json.loads(r'''{params.samples_json}''')

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
# 3) QUAST summary
# -------------------------
import re as _re

def _parse_quast_report(path):
    stats = {{}}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("All statistics"):
                continue
            parts = _re.split(r"\s{{2,}}|\t", line.strip())
            if len(parts) >= 2:
                stats[parts[0].strip()] = parts[-1].strip()
    return stats

_QUAST_KEYS = [
    "# contigs",
    "Largest contig",
    "Total length",
    "GC (%)",
    "N50",
    "N75",
    "L50",
    "L75",
    "# N's per 100 kbp",
]

quast_rows = []
for sample in samples:
    report_path = os.path.join(quast_dir, sample, "report.txt")
    row = {{"sample": sample}}
    if os.path.exists(report_path):
        try:
            stats = _parse_quast_report(report_path)
            for k in _QUAST_KEYS:
                row[k] = stats.get(k, "NA")
        except Exception as e:
            row["parse_error"] = str(e)
    quast_rows.append(row)

quast_df = pd.DataFrame(quast_rows)

# -------------------------
# 4) Kraken2 summary
# -------------------------
kraken_summary_csv = os.path.join(kraken_dir, "kraken2_summary.csv")
if os.path.exists(kraken_summary_csv):
    try:
        kraken_summary_df = pd.read_csv(kraken_summary_csv)
    except Exception:
        kraken_summary_df = pd.DataFrame([{{"Sample": "parse_error"}}])
else:
    kraken_summary_df = pd.DataFrame(columns=["Sample", "Top1", "Top2", "Top3", "Unclassified"])

kraken_db_info = os.path.join(kraken_dir, "kraken2_db_info.txt")
kraken_db_rows = []
if os.path.exists(kraken_db_info):
    with open(kraken_db_info) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t", 1)
            if len(parts) == 2:
                kraken_db_rows.append({{"key": parts[0], "value": parts[1]}})
kraken_db_df = pd.DataFrame(kraken_db_rows, columns=["key", "value"])

# -------------------------
# 5) Autocycler metrics table
# -------------------------
if os.path.exists(autocycler_table_tsv):
    try:
        autocycler_table_df = pd.read_csv(autocycler_table_tsv, sep="\t")
        if autocycler_table_df.empty:
            autocycler_table_df = pd.DataFrame(columns=["sample"])
    except Exception:
        autocycler_table_df = pd.DataFrame([{{"name": "parse_error", "source": autocycler_table_tsv}}])
else:
    autocycler_table_df = pd.DataFrame()

# -------------------------
# 6) Per-sample run status
# -------------------------
status_rows = []
for sample in samples:
    status_rows.append({{
        "sample": sample,
        "qc_clean_fastq": os.path.exists(os.path.join(qc_dir, f"{{sample}}.clean.fastq.gz")),
        "autocycler_fasta": os.path.exists(os.path.join(autocycler_dir, sample, "autocycler_out", "consensus_assembly.fasta")),
        "autocycler_gfa": os.path.exists(os.path.join(autocycler_dir, sample, "autocycler_out", "consensus_assembly.gfa")),
        "dnaapler_gfa": os.path.exists(os.path.join(autocycler_dir, sample, "autocycler_out", "consensus_assembly.reoriented.gfa")),
        "dnaapler_fasta": os.path.exists(os.path.join(autocycler_dir, sample, "autocycler_out", "consensus_assembly.reoriented.fasta")),
        "medaka_final_fasta": os.path.exists(os.path.join(medaka_dir, f"{{sample}}_final_assembly.fasta")),
        "quast_report": os.path.exists(os.path.join(quast_dir, sample, "report.txt")),
        "checkm2_report": os.path.exists(os.path.join(checkm2_dir, sample, "quality_report.tsv")),
        "coverage_tsv": os.path.exists(os.path.join(cov_dir, sample, f"{{sample}}_contig_mean_depth.tsv")),
        "kraken2_report": os.path.exists(os.path.join(kraken_dir, f"{{sample}}.kreport")),
        "bakta_done": os.path.exists(os.path.join(bakta_dir, sample, "DONE.txt")),
    }})

status_df = pd.DataFrame(status_rows)

# -------------------------
# 7) Autocycler helper job summary
# -------------------------
job_rows = []
for sample in samples:
    joblog = os.path.join(autocycler_dir, sample, "assemblies", "joblog.tsv")
    if not os.path.exists(joblog):
        job_rows.append({{
            "sample": sample,
            "joblog_found": False,
            "total_jobs": 0,
            "failed_jobs": 0,
            "failed_commands": ""
        }})
        continue

    try:
        dfj = pd.read_csv(joblog, sep="\t")
        if "Exitval" in dfj.columns:
            failed = dfj[dfj["Exitval"] != 0]
            failed_cmds = "; ".join(failed.get("Command", pd.Series(dtype=str)).astype(str).head(5).tolist())
            job_rows.append({{
                "sample": sample,
                "joblog_found": True,
                "total_jobs": int(dfj.shape[0]),
                "failed_jobs": int(failed.shape[0]),
                "failed_commands": failed_cmds
            }})
        else:
            job_rows.append({{
                "sample": sample,
                "joblog_found": True,
                "total_jobs": int(dfj.shape[0]),
                "failed_jobs": None,
                "failed_commands": "Exitval column not found"
            }})
    except Exception as e:
        job_rows.append({{
            "sample": sample,
            "joblog_found": True,
            "total_jobs": None,
            "failed_jobs": None,
            "failed_commands": f"parse_error: {{e}}"
        }})

autocycler_jobs_df = pd.DataFrame(job_rows)

# -------------------------
# 8) Bad FASTQ summary from prefilter
# -------------------------
bad_fastq_tsv = os.path.join(Path(outdir_root).parent, "0_bad_fastqs", "moved_bad_fastqs.tsv")
if os.path.exists(bad_fastq_tsv):
    try:
        bad_fastq_df = pd.read_csv(bad_fastq_tsv, sep="\t")
    except Exception:
        bad_fastq_df = pd.DataFrame([{{"file_name": "parse_error", "failure_reason": bad_fastq_tsv}}])
else:
    bad_fastq_df = pd.DataFrame(columns=["file_name", "failure_reason"])

# -------------------------
# 9) Stage failure log
# -------------------------
fail_log_tsv = os.path.join(outdir_root, "0_failed_samples.tsv")
if os.path.exists(fail_log_tsv):
    try:
        fail_log_df = pd.read_csv(fail_log_tsv, sep="\t")
    except Exception:
        fail_log_df = pd.DataFrame([{{"sample": "parse_error", "stage": "parse_error", "reason": fail_log_tsv}}])
else:
    fail_log_df = pd.DataFrame(columns=["sample", "stage", "reason"])

# -------------------------
# 10) Write Excel
# -------------------------
outdir = os.path.dirname(out_xlsx)
if outdir:
    os.makedirs(outdir, exist_ok=True)

with pd.ExcelWriter(out_xlsx) as writer:
    checkm2_all.to_excel(writer, sheet_name="CheckM2", index=False)
    coverage_summary.to_excel(writer, sheet_name="Coverage", index=False)
    quast_df.to_excel(writer, sheet_name="QUAST", index=False)
    kraken_summary_df.to_excel(writer, sheet_name="Kraken2", index=False)
    kraken_db_df.to_excel(writer, sheet_name="Kraken2DB", index=False)
    autocycler_table_df.to_excel(writer, sheet_name="AutocyclerTable", index=False)
    status_df.to_excel(writer, sheet_name="RunStatus", index=False)
    autocycler_jobs_df.to_excel(writer, sheet_name="AutocyclerJobs", index=False)
    bad_fastq_df.to_excel(writer, sheet_name="BadFASTQ", index=False)
    fail_log_df.to_excel(writer, sheet_name="FailLog", index=False)

print("[OK] wrote:", out_xlsx)
PY
        """
