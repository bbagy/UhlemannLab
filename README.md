# Uhlemann Lab — Core Pipelines

Containerized Snakemake pipelines for sequencing analysis.  
Each pipeline is self-contained with its own `Dockerfile`, `Snakefile`, `run.sh`, and documentation.

---

## Pipelines

| Pipeline | Type | Description | Docs |
|---------|------|-------------|------|
| **PFsnake** | Illumina WGS | *Plasmodium falciparum* WGS pipeline (QC → mapping → variant calling → CNV → drug resistance → report) | [Open](PFsnake/) |
| **shortWGS** | Illumina WGS | Bacterial WGS pipeline for Illumina short reads | [Open](shortWGS/) |
| **longWGS** | ONT WGS | Bacterial WGS pipeline for ONT long reads | [Open](longWGS/) |
| **RNake** | RNA-seq | RNA-seq pipeline (experimental / under development) | [Open](RNake/) |

---

## General philosophy

- Pipelines are designed for reproducibility and shared workstation usage.
- No tool installation is required on the host machine.
- Each pipeline runs inside Docker.
- Each pipeline folder includes a `run.sh` wrapper so users do not need to type long `docker run` commands.

---

## Quick start

1. Enter a pipeline folder  
2. Build the docker image  
3. Run the pipeline using `run.sh`

Example:

```bash
cd PFsnake
docker build -t pf-snake:1.0 .
./run.sh --help
```

---

## Notes

- Reference databases are **not included** in this repository.
- Output directories and sample data should remain outside the repo.

---

## Maintainer

Developed and maintained by **Heekuk Park**.
