longWGS
=======

Focus
-----

ONT bacterial WGS workflow for assembly, polishing, QC, coverage, and annotation.

Entrypoint
----------

``longWGS/Go_longWGS_V1_1.sh``

Core features
-------------

- FASTQ prefilter (`gzip -t`, bad file quarantine, cached checks)
- Snakemake run with lock auto-retry
- Progress monitor
- Optional Bandage image post-processing

Build
-----

.. code-block:: bash

   cd longWGS
   docker build -t longwgs .

Run
---

.. code-block:: bash

   ./longWGS/Go_longWGS_V1_1.sh -i IN -o OUT -d DB -K

