MAGs
====

Status
------

In testing.

Focus
-----

Modular metagenome QC, assembly/binning, and annotation workflows.

Entrypoints
-----------

``Go_MAGs_QC.sh``, ``Go_MAGs_Assembly.sh``, and ``Go_MAGs_Annotation.sh`` on workstations after ``Go_toWorkstation.sh MAGs``.

Workstation layout
------------------

Wrapper scripts and Snakefiles are placed directly under ``heekuk_path``:
``Go_MAGs_QC.smk``, ``Go_MAGs_Assembly.smk``, and ``Go_MAGs_Annotation.smk``.
Dockerfiles are placed under ``heekuk_path/docker/MAGs``.
Helper scripts are placed under ``heekuk_path/scripts``.

Modules
-------

- QC: paired FASTQ QC and host read depletion
- Assembly: MEGAHIT assembly, binning, coverage, and MAG summary
- Annotation: GTDB-Tk, EggNOG, and KofamScan annotation

Key runtime options
-------------------

- ``-i`` input directory
- ``-o`` output directory
- ``-c`` cores
- ``-j`` Snakemake jobs
- ``-n`` dry-run
- ``-K`` keep-going
