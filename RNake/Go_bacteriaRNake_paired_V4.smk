# 20260501 (Heekuk Park)
# delete annotation.txt. only using gff

import pandas as pd
import glob
import os
import re

# Global variables
master_dir = config["project"] + "_RNAseq_output"
READ_DIR = config["read_dir"]


# === 완전 유연한 R1/R2/SAMPLES 탐색 ===
READS_R1 = sorted(glob.glob(os.path.join(READ_DIR, '*_R1_*.f*q.gz')))
if not READS_R1:
    READS_R1 = sorted([
        f for f in glob.glob(os.path.join(READ_DIR, '*.f*q.gz'))
        if '_R2' not in os.path.basename(f)
    ])

SAMPLES = []
SAMPLE_INFO = {}

for r1 in READS_R1:
    sample_name = re.sub(r'_R1.*\.f(ast)?q\.gz$', '', os.path.basename(r1))

    if sample_name in SAMPLES:
        continue

    # 실제 존재하는 R2 찾기
    possible_r2 = [
        r1.replace('_R1_', '_R2_'),
        r1.replace('_R1.', '_R2.'),
        re.sub(r'_R1([^/]*)\.fastq\.gz$', '_R2\\1.fastq.gz', r1),
        re.sub(r'_R1([^/]*)\.fq\.gz$', '_R2\\1.fq.gz', r1),
    ]

    r2 = next((f for f in possible_r2 if os.path.exists(f)), None)

    if r2:
        SAMPLE_INFO[sample_name] = 'paired'
    else:
        SAMPLE_INFO[sample_name] = 'single'

    SAMPLES.append(sample_name)

# 출력 디렉토리 자동 생성
os.makedirs(master_dir, exist_ok=True)
for subdir in ["1_trim", "2_bowtie2_index", "3_bowtie2_files", "4_htseq-count"]:
    os.makedirs(os.path.join(master_dir, subdir), exist_ok=True)

def detect_feature_type(gff_path, max_scan=200000):
    counts = {}

    with open(gff_path) as f:
        for i, line in enumerate(f):
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            ftype = parts[2]
            counts[ftype] = counts.get(ftype, 0) + 1

            if i >= max_scan:
                break

    if not counts:
        raise ValueError(f"No valid feature types found in {gff_path}")
        
    if "CDS" in counts:
        return "CDS"
    elif "gene" in counts:
        return "gene"
    else:
        return max(counts, key=counts.get)

def detect_best_idattr(gff_path, feature_type=None, max_scan=200000):
    """
    Detect best attribute for htseq-count --idattr.
    Priority: locus_tag > ID > gene > Name > Parent
    """

    if feature_type is None:
        feature_type = detect_feature_type(gff_path)

    priority = ["locus_tag", "ID", "gene", "Name", "Parent"]
    counts = {k: 0 for k in priority}
    total = 0

    with open(gff_path) as f:
        for i, line in enumerate(f):
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")

            # 🔥 FIX
            if len(parts) < 9:
                continue

            if parts[2] != feature_type:
                continue

            total += 1
            attr = parts[8]

            for key in priority:
                if re.search(rf"(?:^|;){key}=", attr):
                    counts[key] += 1

            if total >= 5000 or i >= max_scan:
                break

    if total == 0:
        raise ValueError(f"No feature type '{feature_type}' found in {gff_path}")

    for key in priority:
        if counts[key] / total >= 0.80:
            return key

    return max(priority, key=lambda k: counts[k])


rule all:
    input:
        expand(f"{master_dir}/4_htseq-count/{{sample}}.gene_id.minqual8.txt", sample=SAMPLES),
        f"{master_dir}/2_bowtie2_index/index_build.done",
        f"{master_dir}/merged_counts.csv",
        f"{master_dir}/merged_counts_with_gene_names.csv"

# === helper function ===
def get_r1(wildcards):
    files = glob.glob(os.path.join(READ_DIR, f"{wildcards.sample}_R1*.fastq.gz"))
    if not files:
        # R1 패턴이 없으면 일반 fastq 탐색
        files = glob.glob(os.path.join(READ_DIR, f"{wildcards.sample}*.fastq.gz"))
    if not files:
        raise ValueError(f" R1 FASTQ not found for {wildcards.sample}")
    return files[0]
    
def get_r2(wildcards):
    files = glob.glob(os.path.join(READ_DIR, f"{wildcards.sample}_R2*.fastq.gz"))
    return files if files else []


# === rule trim_reads ===
rule trim_reads:
    input:
        r1=lambda wc: get_r1(wc),
        r2=lambda wc: get_r2(wc)
    params:
        is_paired=lambda wc: len(get_r2(wc)) > 0
    output:
        trimmed_r1_paired = f"{master_dir}/1_trim/{{sample}}.R1.paired.output.fastq.gz",
        trimmed_r2_paired = f"{master_dir}/1_trim/{{sample}}.R2.paired.output.fastq.gz",
        trimmed_r1_unpaired = f"{master_dir}/1_trim/{{sample}}.R1.unpaired.output.fastq.gz",
        trimmed_r2_unpaired = f"{master_dir}/1_trim/{{sample}}.R2.unpaired.output.fastq.gz"
    log:
        f"{master_dir}/1_trim/{{sample}}.trimmomatic.log"
    shell:
        '''
        echo "Trimming {wildcards.sample} ..."

        if [ "{params.is_paired}" = "True" ]; then
            echo "Detected paired-end"
            trimmomatic PE -threads 10 -phred33 \
              {input.r1} {input.r2} \
              {output.trimmed_r1_paired} {output.trimmed_r1_unpaired} \
              {output.trimmed_r2_paired} {output.trimmed_r2_unpaired} \
              LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 2> {log}

        else
            echo "Detected single-end"
            trimmomatic SE -threads 10 -phred33 \
              {input.r1} \
              {output.trimmed_r1_paired} \
              LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 2> {log}

            # dummy outputs
            touch {output.trimmed_r2_paired} {output.trimmed_r2_unpaired} {output.trimmed_r1_unpaired}
        fi
        '''


rule index_genome:
    output:
        touch(f"{master_dir}/2_bowtie2_index/index_build.done")
    params:
        genome=config["genome"]
    shell:
        '''
        bowtie2-build --threads 6 -f {params.genome} {master_dir}/2_bowtie2_index/index
        touch {output}  # Create a dummy file after bowtie2-build is done.
        '''

rule map_reads:
    input:
        r1 = f"{master_dir}/1_trim/{{sample}}.R1.paired.output.fastq.gz",
        r2 = f"{master_dir}/1_trim/{{sample}}.R2.paired.output.fastq.gz",
        index = f"{master_dir}/2_bowtie2_index/index_build.done"
    params:
        is_paired=lambda wc: len(get_r2(wc)) > 0
    output:
        f"{master_dir}/3_bowtie2_files/{{sample}}.aligned.sam"
    shell:
        '''
        echo "Mapping {wildcards.sample} ..."

        if [ "{params.is_paired}" = "True" ]; then
            echo "Detected paired-end"
            bowtie2 --no-unal -p 12 -x {master_dir}/2_bowtie2_index/index \
              -1 {input.r1} -2 {input.r2} \
              -S {output} 2>{master_dir}/3_bowtie2_files/{wildcards.sample}.log
        else
            echo "Detected single-end"
            bowtie2 --no-unal -p 12 -x {master_dir}/2_bowtie2_index/index \
              -U {input.r1} \
              -S {output} 2>{master_dir}/3_bowtie2_files/{wildcards.sample}.log
        fi
        '''

rule sort_sam:
    input:
        f"{master_dir}/3_bowtie2_files/{{sample}}.aligned.sam"
    output:
        f"{master_dir}/3_bowtie2_files/{{sample}}.sorted.sam"
    shell:
        '''
        echo "Processing {wildcards.sample} ..."
        samtools sort -n {input} -o {output}
        '''

rule count_htseq:
    input:
        sam=f"{master_dir}/3_bowtie2_files/{{sample}}.sorted.sam"
    output:
        f"{master_dir}/4_htseq-count/{{sample}}.gene_id.minqual8.txt"
    params:
        reference=config["gff"],
        feature_type=lambda wildcards: detect_feature_type(config["gff"]),
        idattr="ID"
    shell:
        '''
        echo "Processing {wildcards.sample} ..."
        echo "Using feature_type={params.feature_type}, idattr={params.idattr}"

        htseq-count --order=name --stranded=no \
          --type={params.feature_type} \
          --idattr={params.idattr} -a 8 \
          -o {master_dir}/4_htseq-count/{wildcards.sample}.htseq.sam \
          {input.sam} {params.reference} > {output}
        '''


rule merge_counts:
    input:
        files=expand(f"{master_dir}/4_htseq-count/{{sample}}.gene_id.minqual8.txt", sample=SAMPLES)
    output:
        f"{master_dir}/merged_counts.csv"
    run:
        print(f"Merging counts. Input files: {input.files}, Output file: {output}")

        # Create a dictionary where each key is a filename (without extension) and each value is a Series of counts
        dfs = {}
        for file in input.files:
            sample_name = os.path.basename(file).replace(".gene_id.minqual8.txt", "")
            df = pd.read_csv(file, sep='\\t', index_col=0, header=None, names=[sample_name])
            dfs[sample_name] = df[sample_name]


        # Concatenate all Series along the column axis into a DataFrame
        merged_df = pd.concat(dfs, axis=1)

        # Fill any missing values with 0
        merged_df.fillna(0, inplace=True)

        # Write the DataFrame to a new CSV file
        merged_df.to_csv(output[0])

rule map_gene_ids:
    input:
        count_file = master_dir + "/merged_counts.csv",
        annotation_file = config["gff"]
    output:
        mapped_file = master_dir + "/merged_counts_with_gene_names.csv"
    shell:
        r'''
        python - <<'PY'
import pandas as pd
import re

count_df = pd.read_csv("{input.count_file}", index_col=0)
gff = pd.read_csv("{input.annotation_file}", sep="\t", comment="#", header=None)

if (gff[2] == "CDS").sum() > 0:
    gff = gff[gff[2] == "CDS"].copy()
else:
    gff = gff[gff[2] == "gene"].copy()

gff["ID"] = gff[8].str.extract(r"(?:^|;)ID=([^;]+)")
gff["locus_tag"] = gff[8].str.extract(r"(?:^|;)locus_tag=([^;]+)")
gff["gene"] = gff[8].str.extract(r"(?:^|;)gene=([^;]+)")
gff["Name"] = gff[8].str.extract(r"(?:^|;)Name=([^;]+)")
gff["product"] = gff[8].str.extract(r"(?:^|;)product=([^;]+)")

ids = count_df.index.astype(str)

candidates = ["locus_tag", "ID", "gene", "Name"]
best = None
best_hits = -1

for c in candidates:
    if gff[c].isna().all():
        continue
    s = set(gff[c].dropna().astype(str))
    hits = sum(x in s for x in ids)
    if hits > best_hits:
        best_hits = hits
        best = c

if best is None:
    count_df.reset_index(inplace=True)
    count_df.columns.values[0] = "feature_id"
    count_df.insert(1, "symbol", count_df["feature_id"])
    count_df.to_csv("{output.mapped_file}", index=False)
    raise SystemExit

gff = gff.dropna(subset=[best])
gff["symbol"] = gff["gene"].fillna(gff["Name"]).fillna(gff["product"]).fillna(gff[best])

mapping = gff[[best, "symbol"]].drop_duplicates()
mapping = mapping.drop_duplicates(subset=[best], keep="first")

count_df.reset_index(inplace=True)
count_df.columns.values[0] = best
count_df["symbol"] = count_df[best].map(mapping.set_index(best)["symbol"])
count_df["symbol"] = count_df["symbol"].fillna(count_df[best])

cols = count_df.columns.tolist()
rearranged = [cols[0], "symbol"] + [c for c in cols[1:] if c != "symbol"]
count_df = count_df[rearranged]

count_df.to_csv("{output.mapped_file}", index=False)
PY
        '''
