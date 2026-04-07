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
- `KBracken`: Kraken2 + Bracken profiling and merged MPA-style tables
- `MAGs`: modular metagenome QC, MAG assembly/binning, and annotation
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
- Workstation updates are sent with `myscripts/Go_shotgun/Go_toWorkstation.sh`.
- Remote workstation wrappers and Snakefiles are stored under `heekuk_path` with version-free names, while Dockerfiles are stored under `heekuk_path/docker/<pipeline>/`.
- Common wrapper flags:
  - `-n`: dry-run (show execution plan only)
  - `-K`: keep-going (continue independent jobs even if some fail)

---

## Workstation Update Helper

```bash
Go_toWorkstation.sh longWGS
Go_toWorkstation.sh shortWGS
Go_toWorkstation.sh KBracken
Go_toWorkstation.sh MAGs
Go_toWorkstation.sh all
```

The helper sends the latest local versioned workflow files, but stores them on each workstation with stable names:

```text
heekuk_path/Go_longWGS.smk
heekuk_path/Go_shortWGS.smk
heekuk_path/Go_KBracken.smk
heekuk_path/Go_MAGs_QC.smk
heekuk_path/Go_MAGs_Assembly.smk
heekuk_path/Go_MAGs_Annotation.smk
```

---

## Notes

- Reference databases are not versioned in this repository.
- Large outputs should stay outside Git-tracked paths.

---

## Maintainer

Heekuk Park
