# syntax=docker/dockerfile:1
FROM mambaorg/micromamba:1.5.10

ENV MAMBA_DOCKERFILE_ACTIVATE=1 \
    DEBIAN_FRONTEND=noninteractive \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

# -----------------------------
# system packages
# -----------------------------
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
      bash coreutils findutils gawk sed grep \
      git rsync pigz procps \
      parallel \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

# -----------------------------
# Env 1: main pipeline (python 3.11)
# -----------------------------
RUN micromamba create -y -n longwgs -c conda-forge -c bioconda \
      python=3.11 \
      snakemake-minimal=8.* \
      pandas openpyxl biopython \
      samtools minimap2 kraken2 \
      racon raven-assembler \
      porechop_abi nanofilt filtlong dnaapler \
      medaka=2.1.1 \
      flye miniasm minipolish plassembler autocycler \
      bakta=1.12.0 blast=2.17.* \
    && micromamba clean -a -y

# Fail the image build early if key helper tools are not available.
RUN micromamba run -n longwgs autocycler --help >/dev/null \
 && micromamba run -n longwgs minipolish --help >/dev/null \
 && micromamba run -n longwgs filtlong --help >/dev/null \
 && micromamba run -n longwgs dnaapler --help >/dev/null
 

# -----------------------------
# Env 2: checkm2 (python 3.12)
# -----------------------------
RUN micromamba create -y -n checkm2_110 -c conda-forge -c bioconda \
      python=3.12 \
      checkm2=1.1.0 \
    && micromamba clean -a -y

# -----------------------------
# Env 3: QUAST (python 3.10)
#   - quast cannot live with python 3.11
# -----------------------------
RUN micromamba create -y -n quast_env -c conda-forge -c bioconda \
      python=3.10 \
      quast=5.2.0 \
      blast=2.15.* \
    && micromamba clean -a -y

SHELL ["/bin/bash", "-lc"]

# default PATH = longwgs
ENV PATH=/opt/conda/envs/longwgs/bin:$PATH

# DB hints (optional)
ENV BAKTA_DB=/db/bakta_DB
ENV CHECKM2DB=/db/CheckM2_database/uniref100.KO.1.dmnd

WORKDIR /work
CMD ["/bin/bash"]
