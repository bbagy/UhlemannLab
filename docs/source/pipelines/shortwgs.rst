shortWGS
========

Focus
-----

Illumina bacterial WGS workflow for QC, species inference, typing, and summary reporting.

Entrypoints
-----------

- ``shortWGS/Go_shortWGS.sh``
- ``shortWGS/Go_shortWGS_V9_docker.smk``

Core features
-------------

- Input prefilter for FASTQ integrity and paired-read consistency
- Kraken2 classification
- SRST2-based MLST/ARG/Plasmid typing
- TETyper (patched script) integration
- MultiQC + summary report generation

Docker env layout
-----------------

- ``wgs``: main workflow tools (snakemake, fastp, kraken2, R, etc.)
- ``srst2``: SRST2 command execution
- ``tetyper``: TETyper dependencies and patched script runtime

Build
-----

.. code-block:: bash

   cd shortWGS
   docker build --network=host -t shortwgs:1.0 .

Run
---

.. code-block:: bash

   ./shortWGS/Go_shortWGS.sh -i IN -o OUT -d DB -k KRAKEN -r GOWGS -K

