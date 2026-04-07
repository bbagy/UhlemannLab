FROM mambaorg/micromamba:1.5.10

LABEL author="Heekuk Park"
LABEL project="MAGs-QC"

ENV ENV_NAME=mags-qc

RUN micromamba create -y -n ${ENV_NAME} \
    -c conda-forge -c bioconda \
    python=3.11 \
    snakemake=8.25.3 \
    fastp \
    bowtie2 \
    samtools \
    pigz \
    coreutils \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/envs/${ENV_NAME}/bin:$PATH
ENV CONDA_DEFAULT_ENV=${ENV_NAME}

WORKDIR /work
CMD ["bash"]
