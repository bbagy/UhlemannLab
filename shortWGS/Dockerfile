FROM mambaorg/micromamba:1.5.8

USER root

# -----------------------------
# 기본 시스템 패키지
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
# micromamba 환경 생성
# -----------------------------
# 하나의 env에 전부 넣고 PATH 기반으로 통합
RUN micromamba create -y -n wgs \
    -c conda-forge -c bioconda \
    python=3.10 \
    snakemake=8.25.3 \
    fastp \
    bowtie2 \
    samtools \
    bcftools \
    kraken2 \
    multiqc \
    srst2 \
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
# PATH 활성화
# -----------------------------
ENV PATH=/opt/conda/envs/wgs/bin:$PATH
ENV CONDA_DEFAULT_ENV=wgs

# -----------------------------
# TETyper 설치
# (pip로 설치하거나 git clone 방식 가능)
# -----------------------------
RUN pip install --no-cache-dir TETyper

# -----------------------------
# 작업 디렉토리
# -----------------------------
WORKDIR /work

# -----------------------------
# 기본 실행: bash
# -----------------------------
CMD ["bash"]