# UhlemannLab Pipelines

![Snakemake](https://img.shields.io/badge/Snakemake-Workflow-039be5)
![Docker](https://img.shields.io/badge/Docker-Containerized-0db7ed)
![Status](https://img.shields.io/badge/Status-Active-2e7d32)

Containerized sequencing pipelines for routine analysis in the Uhlemann Lab.
Each pipeline is self-contained with its own `Dockerfile`, workflow, wrapper script, and README.

---

## Pipeline Catalog

| Pipeline | Data type | Focus | Entry point |
|---|---|---|---|
| [`PFsnake`](PFsnake/) | Illumina WGS | *P. falciparum* variant/CNV/drug summary | `PFsnake/Go_PFsnake.sh` |
| [`shortWGS`](shortWGS/) | Illumina WGS | Bacterial typing (MLST/ARG/Plasmid/TETyper) | `shortWGS/Go_shortWGS.sh` |
| [`longWGS`](longWGS/) | ONT WGS | Assembly/polishing/QC/annotation | `longWGS/Go_longWGS.sh` |
| [`RNake`](RNake/) | RNA-seq | Bacterial RNA-seq count workflow | `RNake/Go_Rnake.sh` |
| [`TnSeq_ONT`](TnSeq_ONT/) | ONT amplicon/TnSeq | Anchor-based flank extraction + clustering | `TnSeq_ONT/Go_search_tnseq_flank.py` |

---

## Toolset Page

- HTML overview page: [`docs/index.html`](docs/index.html)
- For GitHub Pages, publish from `main` branch / `docs` folder.
- For Read the Docs, this repository now includes:
  - [`.readthedocs.yaml`](.readthedocs.yaml)
  - [`docs/requirements.txt`](docs/requirements.txt)
  - [`docs/source/index.rst`](docs/source/index.rst)

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
