KBracken
========

Status
------

In testing.

Focus
-----

Kraken2 + Bracken profiling with merged MPA-style outputs and master logs.

Entrypoint
----------

``Go_KBracken.sh`` on workstations after ``Go_toWorkstation.sh KBracken``.

Workstation layout
------------------

``Go_KBracken.sh`` and ``Go_KBracken.smk`` are placed directly under ``heekuk_path``.
The Dockerfile is placed under ``heekuk_path/docker/KBracken/Dockerfile``.
Helper scripts are placed under ``heekuk_path/scripts``.

Input naming
------------

Paired-end inference requires explicit ``R1/R2`` or ``forward/reverse`` tokens.
For example, ``ASB_1_R1.fastq.gz`` and ``ASB_1_R2.fastq.gz`` become sample ``ASB_1``.
Files such as ``S1_1.fastq.gz`` and ``S1_2.fastq.gz`` remain separate single-end samples.

Key runtime options
-------------------

- ``-i`` FASTQ directory
- ``-o`` output directory
- ``-d`` Kraken2 database
- ``-c`` cores
- ``-j`` Snakemake jobs
- ``-n`` dry-run
- ``-K`` keep-going
