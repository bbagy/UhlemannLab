FROM mambaorg/micromamba:1.5.8

USER root

# -----------------------------
# Base system packages
# -----------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    wget \
    curl \
    ca-certificates \
    unzip \
    bzip2 \
    procps \
    gawk \
    sed \
    grep \
    coreutils \
    gfortran \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libfontconfig1 \
    libfreetype6 \
    libpng16-16 \
    libharfbuzz0b \
    libfribidi0 \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------
# Env 1: main workflow tools (python 3.11)
# -----------------------------
RUN micromamba create -y -n wgs \
    -c conda-forge -c bioconda \
    python=3.11 \
    snakemake=8.25.3 \
    fastp \
    bowtie2 \
    samtools \
    bcftools \
    kraken2 \
    multiqc \
    blast \
    spades \
    r-base=4.3 \
    r-essentials \
    r-data.table \
    r-tidyverse \
    r-rmarkdown \
    pandoc \
    openjdk \
    && micromamba clean -a -y

# -----------------------------
# Env 2: SRST2 tools + pip-installed SRST2 (python 3.11)
# -----------------------------
RUN micromamba create -y -n srst2 \
    -c conda-forge -c bioconda \
    python=3.11 \
    bowtie2 \
    samtools \
    bcftools \
    && micromamba clean -a -y
RUN micromamba run -n srst2 pip install --no-cache-dir git+https://github.com/katholt/srst2.git

# -----------------------------
# Env 3: TETyper_modi dependencies (python 3.9)
# -----------------------------
RUN micromamba create -y -n tetyper \
    -c conda-forge -c bioconda \
    python=3.9 \
    spades \
    samtools \
    bcftools \
    blast \
    biopython \
    pysam \
    pyvcf3 \
    && micromamba clean -a -y

# -----------------------------
# Default PATH = wgs
# -----------------------------
ENV PATH=/opt/conda/envs/wgs/bin:$PATH
ENV CONDA_DEFAULT_ENV=wgs

# -----------------------------
# Bundle local TETyper patched script
# -----------------------------
COPY TETyper_modi.py /usr/local/bin/TETyper_modi.py
RUN chmod 755 /usr/local/bin/TETyper_modi.py

# -----------------------------
# Working directory
# -----------------------------
WORKDIR /work

# -----------------------------
# Default command: bash
# -----------------------------
CMD ["bash"]
