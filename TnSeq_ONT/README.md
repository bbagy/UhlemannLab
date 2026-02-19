# Tnseq analysis (Anchor → Flank extraction) from ONT reads

---

## Overview

This pipeline extracts **flanking sequences** downstream of a known **anchor sequence**  
(e.g. transposon / cassette / primer junction), using ONT reads.

It supports:
- multiple `.fastq.gz` files in a folder
- optional NanoFilt QC
- optional cutadapt trimming
- minimap2 mapping to the anchor
- flank extraction + junction extraction
- clustering (vsearch)
- summary output (Excel if possible)

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


## Docker run example

```bash
docker run --rm \
  -v /Users/heekukpark/Dropbox/04_scripts/myscripts/Docker/flank_pipeline:/pipeline \
  -v /Users/heekukpark/Documents/03_Core/2026/20260212_Gen_ont:/data \
  -w /pipeline \
  flank-pipeline:1.0 \
  python Go_search_tnseq_flank.py \
    -i /data/fastqs \
    -o /data/output_files \
    --anchor TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCGGGGACTTATCAGCCAACCTGTTAG \
    --flank_len 1000 \
    --threads 8 \
    --cluster_id 0.95
```

---

## Inputs

### `-i / --input`
A directory containing ONT FASTQ files:

- `*.fastq.gz`

Each file is treated as one sample.

Example:
- `KpGD42_barcode01.fastq.gz`
- `KpGD42_barcode02.fastq.gz`

---

## Output directory structure

All outputs are written under the folder specified by `-o`.

Example:

```
output_files/
  anchor.fa
  summary.xlsx
  1_QC/
  2_minimap2/
  3_flank/
  4_cluster/
  logs/
```

---

## Output files (detailed)

### Top-level outputs

- **`anchor.fa`**  
  The anchor sequence written as a FASTA file.  
  (Generated automatically from `--anchor`.)

- **`summary.xlsx`** *(if `openpyxl` is available)*  
  Summary workbook containing:
  - `ReadStats` sheet
  - `TopClusters` sheet

If `openpyxl` is missing, two TSV files are written instead:
- `summary.ReadStats.tsv`
- `summary.TopClusters.tsv`

---

### `1_QC/` (QC + trimming)

For each input FASTQ:

- **`{sample}.clean.fastq.gz`**  
  Cleaned reads after NanoFilt + optional cutadapt trimming.

---

### `2_minimap2/` (anchor mapping)

For each sample:

- **`2_minimap2/{sample}/{sample}.anchor.paf`**  
  minimap2 mapping results in PAF format.

---

### `3_flank/` (flank extraction)

For each sample:

- **`{sample}.flank.fa`**  
  Extracted flank-only sequences:
  - sequence downstream of the anchor hit end
  - anchor itself is NOT included
  - poly-C tail trimming is applied (>=12 C's)
  - sequences shorter than 20 bp are excluded

- **`{sample}.anchorHit_plus_flank.fa`**  
  Junction sequences:
  - includes the anchor-hit region exactly as in the read
  - plus downstream flank region
  - NOT poly-C trimmed
  - sequences shorter than 40 bp are excluded

This file is useful for:
- checking ONT error rate in the anchor region
- confirming exact junction sequence

---

### `4_cluster/` (vsearch clustering)

For each sample:

- **`4_cluster/{sample}/{sample}.centroids.fa`**  
  Representative flank sequences for each cluster.

- **`4_cluster/{sample}/{sample}.clusters.uc`**  
  vsearch UC cluster assignment file.

- **`4_cluster/{sample}/{sample}.cluster_sizes.tsv`**  
  Cluster size table:

  Columns:
  - `cluster_id`
  - `n_reads`

---

### `logs/` (logs + stats)

For each sample:

- **`{sample}.qc_stats.txt`**  
  QC stats table:

  - `raw_reads`
  - `nanofilt_reads`
  - `clean_reads`

- **`{sample}.clean.fastq.flank_stats.txt`**  
  Flank extraction stats:

  - `total_reads`
  - `anchor_hits`
  - `flank_written`
  - `anchorHit_plus_flank_written`
  - `flank_polyC_trimmed`

- **`{sample}.nanofilt.log`** *(only if NanoFilt used)*  
  NanoFilt command log.

- **`{sample}.cutadapt.log`** *(only if cutadapt used)*  
  cutadapt trimming log.

- **`{sample}.clean.fastq.minimap2.log`**  
  minimap2 log.

- **`{sample}.flank.vsearch.log`**  
  vsearch clustering log.

---

## Summary tables explained

### Sheet 1: `ReadStats`

Per sample, includes:

- FASTQ paths
- raw/nanofilt/clean read counts
- QC keep rate
- anchor hit count + hit rate
- number of sequences written to `flank.fa` and `anchorHit_plus_flank.fa`
- length statistics (min/mean/median/max) for both FASTA outputs

---

### Sheet 2: `TopClusters`

Per sample, includes the top N clusters (default: 3):

- cluster_id
- n_reads
- fraction (%)
- centroids fasta path

---

## Notes

- ONT reads may align in reverse orientation; the pipeline automatically reverse-complements reads as needed.
- By default, the pipeline trims poly-C tails (>=12 C's) from flank-only sequences.
- If you want to overwrite existing results, add `--force`.

---


## Maintainer / Contact

Maintained by **Heekuk Park**  
Email: hp2523@cumc.columbia.edu
