# UhlemannLab Pipelines

![Snakemake](https://img.shields.io/badge/Snakemake-Workflow-039be5)
![Docker](https://img.shields.io/badge/Docker-Containerized-0db7ed)
![Status](https://img.shields.io/badge/Status-Active-2e7d32)

Containerized sequencing pipelines for routine analysis in the Uhlemann Lab.
Each pipeline is self-contained with its own `Dockerfile`, workflow, wrapper script, and README.

---

## Pipeline Catalog

- `longWGS`: ONT assembly/polishing/QC/annotation
- `shortWGS`: Illumina typing (MLST/ARG/Plasmid/TETyper)
- `PFsnake`: *P. falciparum* variant/CNV/drug summary
- `RNake`: bacterial RNA-seq count workflow
- `TnSeq_ONT`: ONT flank extraction and clustering

---

## Documentation

- Open full documentation here:
  - **https://<your-project>.readthedocs.io/en/latest/**

---

## Quick Start

```bash
# 1) choose a pipeline
cd PFsnake

# 2) build image
docker build -t pf-snake:1.0 .

# 3) check plan only (dry-run)
./Go_PFsnake.sh -i /path/to/fastq -o out -d /path/to/ref.fasta -n

# 4) run for real
./Go_PFsnake.sh -i /path/to/fastq -o out -d /path/to/ref.fasta
```

---

## Common Conventions

- Data and DB files are mounted from host paths.
- Outputs are written under user-defined output directories.
- Wrapper scripts auto-retry lock issues with Snakemake `--unlock`.
- Dry-run is supported in wrapper scripts with `-n`.

---

## Notes

- Reference databases are not versioned in this repository.
- Large outputs should stay outside Git-tracked paths.

---

## Maintainer

Heekuk Park
