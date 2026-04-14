Quick Start
===========

Prerequisites
-------------

- Docker
- Linux shell
- Pipeline-specific input data and reference databases

Typical flow
------------

1. Move to a pipeline directory.
2. Build the Docker image.
3. Run a dry-run (`-n`) first.
4. Run production with `-K` (keep-going) when appropriate.
5. For longWGS, use `-M permissive` only when you want the summary workbook even if some late-stage jobs fail.

Example (shortWGS)
------------------

.. code-block:: bash

   cd shortWGS
   docker build --network=host -t shortwgs:1.0 .
   ./Go_shortWGS.sh -i /path/fastq -o out -d /path/WGS_DB2 -k /path/kraken -r /path/GoWGS -n
   ./Go_shortWGS.sh -i /path/fastq -o out -d /path/WGS_DB2 -k /path/kraken -r /path/GoWGS -K
