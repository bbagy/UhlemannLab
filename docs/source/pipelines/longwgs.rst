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
- Wrapper option ``-M strict|permissive`` for failure handling
- In ``permissive`` mode, Bakta failures are logged but do not block summary workbook generation

Build
-----

.. code-block:: bash

   cd longWGS
   docker build -t longwgs .

Run
---

.. code-block:: bash

   ./longWGS/Go_longWGS_V1_1.sh -i IN -o OUT -d DB -M strict -K

Permissive mode
---------------

Use ``-M permissive`` when you want the run to continue to the final Excel summary
even if some Bakta jobs fail.

.. code-block:: bash

   ./longWGS/Go_longWGS_V1_1.sh -i IN -o OUT -d DB -M permissive -K
