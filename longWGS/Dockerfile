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
#   - snakemake 8.x needs python >=3.11
#   - medaka works well here
#   - autocycler + assemblers here
# -----------------------------
RUN micromamba create -y -n longwgs -c conda-forge -c bioconda \
      python=3.11 \
      snakemake-minimal=8.* \
      pandas openpyxl biopython \
      samtools minimap2 \
      racon raven-assembler \
      porechop_abi nanofilt \
      medaka=2.1.1 \
      bakta=1.12.0 \
      flye \
      miniasm \
      plassembler \
      autocycler \
    && micromamba clean -a -y

# -----------------------------
# Env 2: checkm2 (python 3.12)
#   - checkm2 1.1.0 wants python 3.12+
# -----------------------------
RUN micromamba create -y -n checkm2_110 -c conda-forge -c bioconda \
      python=3.12 \
      checkm2=1.1.0 \
    && micromamba clean -a -y

# -----------------------------
# default shell
# -----------------------------
SHELL ["/bin/bash", "-lc"]

# default PATH = longwgs env
ENV PATH=/opt/conda/envs/longwgs/bin:$PATH

# DB hints (optional, harmless)
ENV BAKTA_DB=/db/bakta_DB
ENV CHECKM2DB=/db/CheckM2_database/uniref100.KO.1.dmnd

WORKDIR /work
CMD ["/bin/bash"]