# UhlemannLab Pipelines

![Snakemake](https://img.shields.io/badge/Snakemake-Workflow-039be5)
![Docker](https://img.shields.io/badge/Docker-Containerized-0db7ed)
![Status](https://img.shields.io/badge/Status-Active-2e7d32)

Containerized sequencing pipelines for routine analysis in the Uhlemann Lab.
Container launchers share `common/container.sh`; each pipeline keeps its own Dockerfile and workflow.

---

## Pipeline Catalog

Available now (production-ready):

- `longWGS`: ONT assembly/polishing/QC/annotation
- `shortWGS`: Illumina typing (MLST/ARG/Plasmid/TETyper)
- `KBracken`: Kraken2 + Bracken profiling and merged MPA-style tables
- `Humann`: HUMAnN3 + MetaPhlAn4 functional profiling with KEGG orthology output
- `RNake`: bacterial RNA-seq trimming/mapping/count workflow
- `daDake2`: DADA2 amplicon pipeline for 16S (V3V4, V1V2) and ITS with per-project checkpointing

In testing:

- `MAGs`: modular metagenome QC, MAG assembly/binning, and annotation
- `PFsnake`: *P. falciparum* variant/CNV/drug summary
- `TnSeq_ONT`: ONT flank extraction and clustering

---

## Documentation

- Open full documentation here:
  - [Documentation Portal](https://bbagy.github.io/UhlemannLab/)
  - [Docker & Apptainer guide](https://bbagy.github.io/UhlemannLab/containers.html)
  - (Read the Docs URL will be added after first successful RTD build)

The portal serves the static HTML in `docs/`. GitHub Pages must be configured
with **Deploy from a branch → main → /docs**. The `docs/.nojekyll` file keeps
the existing HTML/CSS unchanged during publication. If the GitHub repository
is recreated, restore this Pages setting; pushing the Git history does not
restore repository settings.

---

## Docker / Apptainer selection

For Docker image locations, workstation export, HPC transfer, and SIF conversion,
see the [step-by-step image guide](common/README.md). The reusable helper is
`common/Go_container_image.sh` (`export` on the workstation, `build` on the HPC).

Docker remains the default. The nine container launchers (longWGS, shortWGS,
KBracken, Humann, RNake, PFsnake, and the three MAGs stages) accept:

```bash
# Workstation: existing commands still work.
./Go_shortWGS.sh --container docker [existing options]

# HPC: use a prepared SIF built from the same Docker image.
./Go_shortWGS.sh --container apptainer \
  --container-image /shared/containers/shortwgs.sif [existing options]
```

`--container-image` overrides the image selection; existing `-m` image options
also accept a SIF path for Apptainer. longWGS uses `--container-image` because
it has no `-m` image option. Both `--container=value` and separate-value forms
are supported. No automatic runtime switching or image downloading occurs.

Prepare the SIF once on a machine where pulling/building is permitted:

```bash
# Replace registry/image:tag with the published version of your Docker image.
apptainer pull shortwgs.sif docker://registry/image:tag
```

Local Docker images must first be published or exported for conversion; an
HPC cannot access the workstation's Docker daemon or image cache. Use the
same image version and CPU architecture as the target HPC. Dockerfiles and
Snakefiles are shared; Apptainer translates mounts to bind mounts, runs as the
calling user, and receives explicit environment variables with `--cleanenv`.

When copying launchers manually, copy `common/container.sh` too. In the source
tree it is at `../common/container.sh` relative to each launcher; for flattened
`heekuk_path` installations, place it at `heekuk_path/common/container.sh`.
The workstation update helper must transfer this shared file before launchers.

This adds container execution, not Slurm job distribution. Run inside an HPC
allocation with appropriate cores/memory; rule-level scheduling requires a
separate Snakemake executor/profile setup.

daDake2 and file merge/rename utilities do not currently use containers and
keep their existing host execution. TnSeq_ONT has a Python entry point rather
than a shell launcher; its README includes the direct Apptainer command.

---

## Common Conventions

- Data and DB files are mounted from host paths.
- Outputs are written under user-defined output directories.
- Wrapper scripts auto-retry lock issues with Snakemake `--unlock`.
- Workstation updates are sent with `myscripts/Go_shotgun/Go_toWorkstation.sh`.
- Remote workstation wrappers and Snakefiles are stored under `heekuk_path` with version-free names, while Dockerfiles are stored under `heekuk_path/docker/<pipeline>/`.
- Common wrapper flags:
  - `-n` (or `-x` for `Go_Humann.sh`): dry-run (show execution plan only)
  - `-K`: keep-going (continue independent jobs even if some fail)

---

## Workstation Update Helper

```bash
Go_toWorkstation.sh longWGS
Go_toWorkstation.sh shortWGS
Go_toWorkstation.sh KBracken
Go_toWorkstation.sh Humann
Go_toWorkstation.sh RNake
Go_toWorkstation.sh daDake2
Go_toWorkstation.sh MAGs
Go_toWorkstation.sh all
```

The helper sends the latest local versioned workflow files, but stores them on each workstation with stable names:

```text
heekuk_path/Go_longWGS.smk
heekuk_path/Go_shortWGS.smk
heekuk_path/Go_KBracken.smk
heekuk_path/Go_Humann.smk
heekuk_path/Go_RNake.smk
heekuk_path/daDake2/Go_daDake2.smk       ← symlink: heekuk_path/Go_daDake2.sh
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
