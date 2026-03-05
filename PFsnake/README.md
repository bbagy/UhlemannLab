# PFsnake

![Illumina](https://img.shields.io/badge/Reads-Illumina-1f77b4)
![Workflow](https://img.shields.io/badge/Workflow-Snakemake-039be5)
![Container](https://img.shields.io/badge/Runtime-Docker-0db7ed)

*Plasmodium falciparum* WGS workflow for QC, mapping, variant calling, CNV, mixed infection metrics, and report generation.

---

## Pipeline

```mermaid
flowchart LR
  A[FASTQ] --> B[fastp + MultiQC]
  B --> C[BWA + Picard]
  C --> D[GATK variants]
  D --> E[CNV / Mixed / Recomb]
  E --> F[Drug summary]
  F --> G[Excel reports]
```

---

## Requirements

- Docker
- Linux environment
- Illumina FASTQ (SE or PE)
- Reference FASTA (`-d`, required)
- PlasmoDB GFF (`-g`, optional but required for gene-CNV/drug-gene outputs)

---

## Build

```bash
docker build -t pf-snake:1.0 .
```

---

## Quick Start

```bash
Go_PFsnake.sh \
  -i /path/to/fastq \
  -o pfsnake_out \
  -d /path/to/Pf3D7.fasta \
  -g /path/to/PlasmoDB.gff \
  -c 8 -t 4 -p 2
```

---

## Options

| Flag | Default | Description |
|---|---:|---|
| `-i` | - | Input FASTQ directory |
| `-o` | - | Output directory |
| `-d` | - | Reference FASTA |
| `-s` | script directory | Optional Snakefile directory override |
| `-g` | (none) | GFF for gene-CNV/drug gene summary |
| `-c` | `8` | Snakemake cores |
| `-t` | `4` | Thread setting for selected rules |
| `-p` | `2` | Read mode (`2` PE, `1` SE) |
| `-m` | `pf-snake:1.0` | Docker image |
| `-n` | off | Dry-run (`snakemake --dry-run`) |
| `-K` | off | Keep going on independent job failures (`--keep-going`) |
| `-P` | `1` | Show periodic progress snapshots in current terminal (`1` on, `0` off) |

---

## Inputs

Supported sample naming:

- `<sample>_R1_001.fastq.gz` / `<sample>_R2_001.fastq.gz`
- `<sample>_R1.fastq.gz` / `<sample>_R2.fastq.gz`
- `<sample>.R1.fastq.gz` / `<sample>.R2.fastq.gz`

---

## Outputs

```text
OUTDIR/
  1_QC/
    fastp/
    multiqc/
  2_bwa/
  3_bam_qc/
  4_gatk_out/
  5_variant_summary/
  6_cnv/
    window/
    gene/                (if -g provided)
  7_mixed_infection/
  8_recomb_explore/
  9_drug_resistance/     (if -g provided)
  logs/
  benchmarks/
  YYYYMMDD_PF_QC_report.xlsx
  YYYYMMDD_PF_variant_summary.xlsx
```

---

## Notes

- `Go_PFsnake.sh` runs input prefilter before Snakemake:
  - gzip integrity check (`gzip -t`)
  - bad FASTQ moved to sibling `0_bad_fastqs/`
  - bad log: `0_bad_fastqs/moved_bad_fastqs.tsv`
  - when `-p 2`, unmatched R1/R2 files are also moved to bad folder

---

## Maintainer

Heekuk Park
