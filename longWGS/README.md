# longWGS — Bacterial WGS (ONT)

---

## Overview

This Snakemake-based pipeline performs bacterial whole-genome assembly and annotation from Oxford Nanopore (ONT) reads.  
It runs QC → assembly (autocycler) → polishing (medaka) → QC (QUAST/CheckM2) → coverage → annotation (Bakta), inside a Docker container.

---

## Requirements

- Docker  
- A Linux server
- Long-read FASTQ files (`*.fastq.gz`)

---

## Build docker image

```bash
docker build -t longwgs:1.0 .
```

---

## Run

```bash
Go_longWGS.sh -i [input DIR] -o [output DIR] -d [DB location] -s [snakemake file location] -p 0

Go_longWGS.sh -i 1_merged_fastqs -o 2_longWGS_out -d /media/uhlemann/Core3_V2/DB/longWGS_DB -s /home/uhlemann/heekuk_path -p 0
```

### Options

- `-p 0` : skip porechop (recommended if your reads are already demultiplexed/trimmed)
- `-p 1` : run porechop_abi before NanoFilt

---

## Inputs

### 1) Input FASTQ directory (`-i`)

A directory containing ONT reads as gzip-compressed FASTQ files:

- `*.fastq.gz` or `*.fq.gz`
- one file per sample (e.g., `barcode01.fastq.gz`)

Example:

```
1_merged_fastqs/
  barcode01.fastq.gz
  barcode02.fastq.gz
  barcode03.fastq.gz
```

⚠️ Notes:
- Corrupted `.fastq.gz` will fail QC (NanoFilt).  
  You can test quickly using:
  ```bash
  gzip -t barcode01.fastq.gz
  ```

---

### 2) DB root directory (`-d`)

The DB folder must contain BOTH Bakta DB and CheckM2 DB under a single root:

```
DB/
  bakta_DB/
  CheckM2_database/
    uniref100.KO.1.dmnd
```

Example:

```bash
-d /media/uhlemann/Core3_V2/DB/longWGS_DB
```

---

### 3) Snakemake directory (`-s`)

This is the folder containing the pipeline Snakefile:

```
/home/uhlemann/heekuk_path/
  Go_longWGS.smk
```

Example:

```bash
-s /home/uhlemann/heekuk_path
```

---

## Outputs

The output directory (`-o`) will be created with the following structure:

```
OUTDIR/
  1_QC/             (porechop_abi + NanoFilt output: *.clean.fastq.gz)
  2_quast/          (QUAST reports; from medaka fasta)
  3_autocycler/     ({sample}/autocycler_out/consensus_assembly.{fasta,gfa})
  4_medaka/         ({sample}_final_assembly.fasta)
  5_checkm2/        ({sample}/quality_report.tsv + DONE.txt + raw/)
  6_coverage/       ({sample}/sorted BAM + contig mean depth + DONE.txt)
  7_bakta/          ({sample}/DONE.txt + annotation outputs)
  checkm2_coverage_summary.xlsx
```

---

## Notes

- Snakemake is executed inside Docker.
- The pipeline automatically retries with `--unlock` if a Snakemake lock issue is detected.
- If you want Bandage images, run Bandage separately after the pipeline finishes  
  (Bandage is intentionally NOT included in the Snakefile).

---

## Maintainer

Heekuk Park
