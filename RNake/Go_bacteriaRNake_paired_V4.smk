# 20230922 (Heekuk Park)
# delete annotation.txt. only using gff

import pandas as pd
import glob
import os
import re

# Global variables
master_dir = config["project"] + "_RNAseq_output"
READ_DIR = config["read_dir"] 
conda_shotgun = "/media/uhlemann/core4/DB/MAGs/shotgun.yaml"
conda_MAGs = "/media/uhlemann/core4/DB/MAGs/MAGs.yaml"


# === 완전 유연한 R1/R2/SAMPLES 탐색 ===
READS_R1 = sorted(glob.glob(os.path.join(READ_DIR, '*R1*.f*q.gz')))
READS_R2 = []
SAMPLES = []

for r1 in READS_R1:
    # 가능한 모든 R2 후보 생성 (_R1_001, _R1_nohuman, _R1trim 등 모든 경우 대응)
    r2_candidates = [
        re.sub(r'_R1([^/]*)\.fastq\.gz$', '_R2\\1.fastq.gz', r1),
        re.sub(r'_R1([^/]*)\.fq\.gz$', '_R2\\1.fq.gz', r1),
        re.sub(r'_R1([^/]*)\.fastq$', '_R2\\1.fastq', r1),
        r1.replace('_R1_', '_R2_'),
        r1.replace('R1', 'R2'),
    ]

    # 존재하는 R2 파일 중 첫 번째 사용
    r2 = next((f for f in r2_candidates if os.path.exists(f)), None)

    if r2:
        READS_R2.append(r2)
        # 샘플 이름 추출 (_R1 앞까지만 자동으로)
        sample_name = re.split(r'_R1[^/]*\.f(ast)?q\.gz$', os.path.basename(r1))[0]
        SAMPLES.append(sample_name)
    else:
        print(f"⚠️ Warning: R2 not found for {r1}")

# 출력 디렉토리 자동 생성
os.makedirs(master_dir, exist_ok=True)
for subdir in ["1_trim", "2_bowtie2_index", "3_bowtie2_files", "4_htseq-count"]:
    os.makedirs(os.path.join(master_dir, subdir), exist_ok=True)


def detect_best_idattr(gff_path, feature_type="CDS", max_scan=200000):
    """
    Detect best attribute for htseq-count --idattr.
    Priority: locus_tag > ID > gene > Name > Parent
    Returns string.
    """
    priority = ["locus_tag", "ID", "gene", "Name", "Parent"]
    counts = {k: 0 for k in priority}
    total = 0

    with open(gff_path) as f:
        for i, line in enumerate(f):
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != feature_type:
                continue

            total += 1
            attr = parts[8]

            for key in priority:
                if re.search(rf"(?:^|;){key}=", attr):
                    counts[key] += 1

            if total >= 5000:
                break
            if i >= max_scan:
                break

    if total == 0:
        raise ValueError(f"No feature type '{feature_type}' found in {gff_path}")

    # choose highest priority that exists in >=80% of CDS
    for key in priority:
        if counts[key] / total >= 0.80:
            return key

    # fallback: choose max coverage
    best = max(priority, key=lambda k: counts[k])
    return best


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
        raise ValueError(f"❌ R1 FASTQ not found for {wildcards.sample}")
    return files[0]

def get_r2(wildcards):
    files = glob.glob(os.path.join(READ_DIR, f"{wildcards.sample}_R2*.fastq.gz"))
    if not files:
        raise ValueError(f"❌ R2 FASTQ not found for {wildcards.sample}")
    return files[0]


# === rule trim_reads ===
rule trim_reads:
    input:
        r1 = lambda wildcards: get_r1(wildcards),
        r2 = lambda wildcards: get_r2(wildcards)
    output:
        trimmed_r1_paired = f"{master_dir}/1_trim/{{sample}}.R1.paired.output.fastq.gz",
        trimmed_r2_paired = f"{master_dir}/1_trim/{{sample}}.R2.paired.output.fastq.gz",
        trimmed_r1_unpaired = f"{master_dir}/1_trim/{{sample}}.R1.unpaired.output.fastq.gz",
        trimmed_r2_unpaired = f"{master_dir}/1_trim/{{sample}}.R2.unpaired.output.fastq.gz"
    log:
        f"{master_dir}/1_trim/{{sample}}.trimmomatic.log"
    conda:
        conda_shotgun
    shell:
        """
        echo "Trimming {wildcards.sample} ..."
        trimmomatic PE -threads 10 -summary {master_dir}/1_trim/{wildcards.sample}.statssummary -phred33 \
          {input.r1} {input.r2} \
          {output.trimmed_r1_paired} {output.trimmed_r1_unpaired} \
          {output.trimmed_r2_paired} {output.trimmed_r2_unpaired} \
          LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 2> {log}
        """


rule index_genome:
    output:
        touch(f"{master_dir}/2_bowtie2_index/index_build.done")
    params:
        genome=config["genome"]
    conda:
        conda_MAGs
    shell:
        """
        bowtie2-build --threads 6 -f {params.genome} {master_dir}/2_bowtie2_index/index
        touch {output}  # Create a dummy file after bowtie2-build is done.
        """



rule map_reads:
    input:
        r1 = f"{master_dir}/1_trim/{{sample}}.R1.paired.output.fastq.gz",
        r2 = f"{master_dir}/1_trim/{{sample}}.R2.paired.output.fastq.gz",
        index = f"{master_dir}/2_bowtie2_index/index_build.done"
    output:
        f"{master_dir}/3_bowtie2_files/{{sample}}.aligned.sam"
    conda:
        conda_MAGs
    shell:
        """
        echo "Mapping {wildcards.sample} against indexed NR5452 WT SPADES assembly..."
        bowtie2 --no-unal -p 12 -x {master_dir}/2_bowtie2_index/index -1 {input.r1} -2 {input.r2} -S {output} 2>{master_dir}/3_bowtie2_files/{wildcards.sample}.log
        """

rule sort_sam:
    input:
        f"{master_dir}/3_bowtie2_files/{{sample}}.aligned.sam"
    output:
        f"{master_dir}/3_bowtie2_files/{{sample}}.sorted.sam"
    conda:
        conda_MAGs
    shell:
        """
        echo "Processing {wildcards.sample} ..."
        samtools sort -n {input} -o {output}
        """

rule count_htseq:
    input:
        sam=f"{master_dir}/3_bowtie2_files/{{sample}}.sorted.sam"
    output:
        f"{master_dir}/4_htseq-count/{{sample}}.gene_id.minqual8.txt"
    params:
        reference=config["gff"],
        idattr=lambda wildcards: detect_best_idattr(config["gff"], "CDS")
    conda:
        conda_shotgun
    shell:
        """
        echo "Processing {wildcards.sample} ..."
        echo "Using idattr={params.idattr}"
        htseq-count --order=name --stranded=no --type=CDS --idattr={params.idattr} -a 8 \
          -o {master_dir}/4_htseq-count/{wildcards.sample}.htseq.sam \
          {input.sam} {params.reference} > {output}
        """


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
            sample_name = file.split('/')[-1].split('.')[0] # Change this line to correctly extract sample name
            df = pd.read_csv(file, sep='\t', index_col=0, header=None, names=[sample_name])
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
    conda:
        conda_shotgun
    shell:
        r"""
        python - <<'PY'
import pandas as pd
import re

count_df = pd.read_csv("{input.count_file}", index_col=0)
gff = pd.read_csv("{input.annotation_file}", sep="\t", comment="#", header=None)

gff = gff[gff[2] == "CDS"].copy()

# Extract all possible IDs
gff["ID"] = gff[8].str.extract(r"(?:^|;)ID=([^;]+)")
gff["locus_tag"] = gff[8].str.extract(r"(?:^|;)locus_tag=([^;]+)")
gff["gene"] = gff[8].str.extract(r"(?:^|;)gene=([^;]+)")
gff["Name"] = gff[8].str.extract(r"(?:^|;)Name=([^;]+)")
gff["product"] = gff[8].str.extract(r"(?:^|;)product=([^;]+)")

# Decide which key matches count_df index best
ids = count_df.index.astype(str)

candidates = ["locus_tag", "ID", "gene", "Name"]
best = None
best_hits = -1

for c in candidates:
    if gff[c].isna().all():
        continue
    s = set(gff[c].dropna().astype(str).tolist())
    hits = sum([1 for x in ids if x in s])
    if hits > best_hits:
        best_hits = hits
        best = c

if best is None:
    # fallback: no mapping possible
    count_df.reset_index(inplace=True)
    count_df.columns.values[0] = "feature_id"
    count_df.insert(1, "symbol", count_df["feature_id"])
    count_df.to_csv("{output.mapped_file}", index=False)
    raise SystemExit

# Make symbol column (priority)
gff = gff.dropna(subset=[best])
gff["symbol"] = gff["gene"].fillna(gff["Name"]).fillna(gff["product"]).fillna(gff[best])

mapping = gff[[best, "symbol"]].drop_duplicates()
mapping = mapping.drop_duplicates(subset=[best], keep="first")

# Apply mapping
count_df.reset_index(inplace=True)
count_df.columns.values[0] = best
count_df["symbol"] = count_df[best].map(mapping.set_index(best)["symbol"])
count_df["symbol"] = count_df["symbol"].fillna(count_df[best])

# Put symbol second
cols = count_df.columns.tolist()
rearranged = [cols[0], "symbol"] + [c for c in cols[1:] if c != "symbol"]
count_df = count_df[rearranged]

count_df.to_csv("{output.mapped_file}", index=False)
PY
        """



