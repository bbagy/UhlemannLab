# MAGs

Modular MAG workflow collection for Uhlemann Lab metagenomics.

The modules are intentionally separated so QC, assembly/binning, inStrain, and annotation can be run independently.

## Layout

```text
MAGs/
  workflow/
    Go_MAGs_QC_V1.smk
    Go_MAGs_Assembly_V1.smk
    Go_MAGs_Annotation_V1.smk
  docker/
    Dockerfile.qc
    Dockerfile.assembly
    Dockerfile.annotation
  scripts/
    summarize_all_results.py
  Go_MAGs_QC.sh
  Go_MAGs_Assembly.sh
  Go_MAGs_Annotation.sh
  MAGs_pipeline.png
  MAGs_pipeline.pdf
```

## Workstation Layout

After `Go_toWorkstation.sh MAGs`, files are placed as:

```text
/home/uhlemann*/heekuk_path/
  Go_MAGs_QC.sh
  Go_MAGs_QC.smk
  Go_MAGs_Assembly.sh
  Go_MAGs_Assembly.smk
  Go_MAGs_Annotation.sh
  Go_MAGs_Annotation.smk
  scripts/
    summarize_all_results.py
  docker/MAGs/
    Dockerfile.qc
    Dockerfile.assembly
    Dockerfile.annotation
```

## QC Module

The QC workflow performs paired-end FASTQ QC and host read depletion:

```text
raw FASTQ
  -> fastp
  -> filtered_fastq/
  -> bowtie2 host filtering
  -> host_filtered_fastq/*_R1_nohuman.fastq.gz
  -> host_filtered_fastq/*_R2_nohuman.fastq.gz
  -> QC_summary.csv
```

### Build QC Image

Build from the workstation Docker directory:

```bash
cd /home/uhlemann/heekuk_path/docker/MAGs
docker build -f Dockerfile.qc -t mags-qc:1.0 .
```

### Build Assembly Image

Build from the workstation Docker directory:

```bash
cd /home/uhlemann/heekuk_path/docker/MAGs
docker build -f Dockerfile.assembly -t mags-assembly:1.0 .
```

### Build Annotation Image

Build from the workstation Docker directory:

```bash
cd /home/uhlemann/heekuk_path/docker/MAGs
docker build -f Dockerfile.annotation -t mags-annotation:1.0 .
```

### Run QC With Docker Wrapper

`-d` must be the Bowtie2 host index prefix, not just the directory.

```bash
./Go_MAGs_QC.sh \
  -i /path/to/input_fastqs \
  -o mags_qc_out \
  -d /path/to/bowtie2_host_index/hg38 \
  -c 8 \
  -j 4 \
  -K
```

### Direct Snakemake Run

```bash
fastq_dir="input_fastqs"
output_dir="mags_qc_out"
host_db="/path/to/bowtie2_host_index/hg38"

snakemake --snakefile /home/uhlemann/heekuk_path/Go_MAGs_QC.smk \
  --config fastq_dir="$fastq_dir" output_dir="$output_dir" host_db="$host_db" paired=2 threads=8 \
  --cores 8 --jobs 4 \
  --latency-wait 60 --rerun-incomplete
```

## QC Input Naming

Supported paired-end R1 naming patterns:

- `<sample>_L001_R1_001.fastq.gz`
- `<sample>_R1_001.fastq.gz`
- `<sample>_R1.fastq.gz`
- `<sample>.R1.fastq.gz`

Multiple lanes are concatenated per sample before `fastp`.

## QC Output

```text
mags_qc_out/
  filtered_fastq/
    <sample>_R1_filtered.fastq.gz
    <sample>_R2_filtered.fastq.gz
    <sample>_fastp.json
    <sample>_fastp.html
  host_filtered_fastq/
    <sample>_R1_nohuman.fastq.gz
    <sample>_R2_nohuman.fastq.gz
  intermediate/
  logs/
  QC_summary.csv
```

## Assembly And Binning Module

After QC is complete, use the QC `host_filtered_fastq` directory as the assembly input:

```bash
./Go_MAGs_Assembly.sh \
  -i mags_qc_out/host_filtered_fastq \
  -o mags_assembly_out \
  -d /path/to/checkm_data \
  -c 24 \
  -j 4 \
  -t 24 \
  -M 128000 \
  -K
```

Dry-run first:

```bash
./Go_MAGs_Assembly.sh \
  -i mags_qc_out/host_filtered_fastq \
  -o mags_assembly_out \
  -d /path/to/checkm_data \
  -c 24 \
  -j 4 \
  -n
```

### Direct Snakemake Run

```bash
snakemake --snakefile /home/uhlemann/heekuk_path/Go_MAGs_Assembly.smk \
  --config fastq_dir="mags_qc_out/host_filtered_fastq" output_dir="mags_assembly_out" \
  checkm_data_dir="/path/to/checkm_data" megahit_threads=24 megahit_memory=128000 binning_tools="concoct,metabat2,maxbin2" \
  --cores 24 --jobs 4 \
  --latency-wait 60 --rerun-incomplete --keep-going
```

### Assembly Output

```text
mags_assembly_out/
  1_MEGAHIT/
  2_Contigs/
  3_Binning/
  4_Mapping/
  5_Coverage/
  6_Mapped_Reads/
  7_DAS_tool_out/
  8_checkM_summary/
  9_Final_MAGs/
  assembly_summary.csv
  checkm_summary_combined.csv
```

## Annotation Module

The annotation workflow is also Docker-based. Provide a directory containing final MAG FASTA files (`*.fa` or `*.fasta`) plus the GTDB-Tk, EggNOG, and KofamScan database directories.

```bash
./Go_MAGs_Annotation.sh \
  -i mags_assembly_out/9_Final_MAGs \
  -o mags_annotation_out \
  -g /path/to/gtdbtk/release220 \
  -e /path/to/eggnog \
  -k /path/to/kofam_scan \
  -c 24 \
  -j 4 \
  -K
```

Dry-run first:

```bash
./Go_MAGs_Annotation.sh \
  -i mags_assembly_out/9_Final_MAGs \
  -o mags_annotation_out \
  -g /path/to/gtdbtk/release220 \
  -e /path/to/eggnog \
  -k /path/to/kofam_scan \
  -c 24 \
  -j 4 \
  -n
```

Expected mounted database layout inside Docker:

```text
/db/gtdbtk/
  mash_db/genomic_mash_db.msh
/db/eggnog/
  eggnog_proteins.dmnd
/db/kofam_scan/
  ko_list
  profiles/
```
