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
  - [Documentation Portal](https://bbagy.github.io/UhlemannLab/)
  - (Read the Docs URL will be added after first successful RTD build)

---

## Common Conventions

- Data and DB files are mounted from host paths.
- Outputs are written under user-defined output directories.
- Wrapper scripts auto-retry lock issues with Snakemake `--unlock`.
- Common wrapper flags:
  - `-n`: dry-run (show execution plan only)
  - `-K`: keep-going (continue independent jobs even if some fail)

---

## Notes

- Reference databases are not versioned in this repository.
- Large outputs should stay outside Git-tracked paths.

---

## Maintainer

Heekuk Park
