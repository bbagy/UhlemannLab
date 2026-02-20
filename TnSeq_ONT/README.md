# TnSeq_ONT

![ONT](https://img.shields.io/badge/Reads-ONT-2ca02c)
![Python](https://img.shields.io/badge/Entry-Python%20CLI-3776ab)
![Container](https://img.shields.io/badge/Runtime-Docker-0db7ed)

Anchor-based flank extraction workflow from ONT reads, with optional QC/trimming, mapping, clustering, and summary tables.

---

## Pipeline

```mermaid
flowchart LR
  A[FASTQ] --> B[Optional NanoFilt/cutadapt]
  B --> C[minimap2 to anchor]
  C --> D[Flank extraction]
  D --> E[vsearch clustering]
  E --> F[summary.xlsx]
```

---

## Requirements

- Docker
- Linux environment
- ONT `*.fastq.gz`

---

## Build

```bash
docker build -t flank-pipeline:1.0 .
```

---

## Quick Start

```bash
docker run --rm \
  -v /path/to/TnSeq_ONT:/pipeline \
  -v /path/to/data:/data \
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

## Key Arguments

| Argument | Default | Description |
|---|---:|---|
| `-i, --input` | - | Input FASTQ directory |
| `-o, --output` | - | Output directory |
| `--anchor` | - | Anchor sequence |
| `--anchor_name` | `anchor` | FASTA header for anchor |
| `--flank_len` | `300` | Extracted downstream length |
| `--threads` | `4` | minimap2 threads |
| `--qmin` | none | NanoFilt min quality |
| `--min_len` | none | NanoFilt min read length |
| `--cluster_id` | `0.95` | vsearch identity |
| `--top_n` | `3` | Top clusters in summary |
| `--force` | off | Overwrite existing outputs |

---

## Outputs

```text
output_files/
  anchor.fa
  summary.xlsx                       (or summary.ReadStats.tsv / summary.TopClusters.tsv)
  1_QC/
  2_minimap2/
  3_flank/
  4_cluster/
  logs/
```

Important per-sample files:

- `3_flank/{sample}.flank.fa`
- `3_flank/{sample}.anchorHit_plus_flank.fa`
- `4_cluster/{sample}/{sample}.cluster_sizes.tsv`
- `logs/{sample}.qc_stats.txt`
- `logs/{sample}.clean.fastq.flank_stats.txt`

---

## Maintainer

Heekuk Park
