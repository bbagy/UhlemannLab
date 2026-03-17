# longWGS

![ONT](https://img.shields.io/badge/Reads-ONT-2ca02c)
![Workflow](https://img.shields.io/badge/Workflow-Snakemake-039be5)
![Container](https://img.shields.io/badge/Runtime-Docker-0db7ed)

Bacterial ONT WGS workflow for assembly, polishing, QC, coverage, and annotation.

## Current Entrypoint

Use `Go_longWGS_V1_1.sh` as the current wrapper.

```bash
./longWGS/Go_longWGS_V1_1.sh -h
```

## Pipeline

```mermaid
flowchart LR
  P[barcode dirs\nfastq_pass/] --> Q[Go_merge_rename.sh\nmerge + rename]
  Q --> A[sample.fastq.gz]
  A --> B[Prefilter gzip check]
  B --> C[QC / optional porechop]
  C --> D[Autocycler assembly]
  D --> E[Medaka polish]
  E --> F[QUAST + CheckM2]
  F --> G[Coverage]
  G --> H[Bakta annotation]
  H --> I[Optional Bandage images]
```

## Requirements

- Docker
- Linux shell
- ONT FASTQ files (`*.fastq.gz` or `*.fq.gz`)
- DB directory mounted as `-d`, including:
  - `medaka_models/`
  - `plassembler_db/`
  - databases used by rules (for example Bakta and CheckM2)

## Build

The wrapper runs image name `longwgs` directly, so tag exactly as below.

```bash
cd longWGS
docker build -t longwgs .
```

## Docker Tool Inventory

Main tools included in the image (pipeline-dependent):

- `snakemake`
- `autocycler`
- `medaka`
- `quast`
- `checkm2`
- `bakta`
- `samtools`
- supporting utilities used by workflow rules

## Run A Single Tool From The Image

Examples (without conda on host):

```bash
# Snakemake version
docker run --rm longwgs snakemake --version

# CheckM2 help (example)
docker run --rm longwgs checkm2 --help

# Bakta help (example)
docker run --rm longwgs bakta --help
```

## Pre-processing: Merge & Rename

ONT sequencing output typically produces per-barcode directories (e.g. `fastq_pass/barcode01/`).
Before running the pipeline, merge the per-read files into a single fastq.gz per sample and
rename them from barcode names to sample names using `Go_merge_rename.sh`.

```text
fastq_pass/
  barcode01/  →  merged_fastqs/KP0011.fastq.gz
  barcode02/  →  merged_fastqs/KP0063.fastq.gz
  ...
```

**Map file format** (tab or space separated, `sample  barcode`):

```text
KP0011  barcode01
KP0063  barcode02
...
```

**Run:**

```bash
# dry-run first
bash Go_merge_rename.sh -i fastq_pass -o merged_fastqs -m samples_barcodes.txt -n

# apply
bash Go_merge_rename.sh -i fastq_pass -o merged_fastqs -m samples_barcodes.txt
```

The `merged_fastqs/` directory is then used as input (`-i`) for the pipeline.

## Quick Start

```bash
./longWGS/Go_longWGS_V1_1.sh \
  -i /path/to/fastq \
  -o /path/to/output \
  -d /path/to/db \
  -p 0 \
  -K
```

Real example:

```bash
Go_longWGS.sh \
  -i 1_merged_fastqs \
  -o 2_longWGS_out \
  -d /media/uhlemann/Core3_V2/DB/longWGS_DB \
  -s /home/uhlemann/heekuk_path \
  -p 0 \
  -K
```

## Options

| Flag | Default | Description |
|---|---:|---|
| `-i` | - | Input FASTQ directory |
| `-o` | - | Output directory |
| `-d` | - | DB root mounted to container as `/db` |
| `-s` | script directory | Optional Snakefile directory override |
| `-p` | `0` | `0`: skip porechop, `1`: run porechop |
| `-n` | off | Dry-run (`snakemake --dry-run`) |
| `-K` | off | Keep going on independent failures (`--keep-going`) |

## What The Wrapper Does

1. Prefilter input FASTQ files before Snakemake.
2. Move corrupted FASTQ files to sibling folder `0_bad_fastqs/`.
3. Cache prefilter state in `0_bad_fastqs/DONE.txt`.
4. Print periodic progress snapshots while running.
5. Auto-retry lock errors (`--unlock` then rerun).
6. Handle `INT/TERM` and clean child/background processes.
7. Generate Bandage images when `Bandage` is available on host and GFA exists.

## Input Layout

```text
FASTQ_DIR/
  sample1.fastq.gz
  sample2.fastq.gz
  ...
```

## Output Layout

```text
OUTDIR/
  1_QC/
  2_quast/
  3_autocycler/
  4_medaka/
  5_checkm2/
  6_coverage/
  7_bakta/
  8_Bandage_image/          # created when Bandage step runs
  0_failed_samples.tsv      # created by workflow on rule failures
```

Prefilter artifacts:

```text
<parent_of_FASTQ_DIR>/0_bad_fastqs/
  DONE.txt
  moved_bad_fastqs.tsv
  *.fastq.gz                # moved bad files
```

## Common Runs

Dry-run:

```bash
./longWGS/Go_longWGS_V1_1.sh -i IN -o OUT -d DB -n
```

Production run (recommended):

```bash
./longWGS/Go_longWGS_V1_1.sh -i IN -o OUT -d DB -K
```

Enable porechop:

```bash
./longWGS/Go_longWGS_V1_1.sh -i IN -o OUT -d DB -p 1 -K
```

## Troubleshooting

- `docker: image not found longwgs`
  - build with `docker build -t longwgs longWGS`
- lock-related failure in log
  - wrapper already retries with `--unlock`
- repeated prefilter not running
  - expected when input files are older than `0_bad_fastqs/DONE.txt`
- no Bandage images
  - install `Bandage` on host or skip (pipeline core output is unaffected)

## Maintainer

Heekuk Park
