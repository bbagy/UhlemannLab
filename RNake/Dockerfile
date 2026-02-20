FROM mambaorg/micromamba:latest

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN micromamba install -y -n base -c conda-forge -c bioconda -c defaults \
    python=3.11 \
    snakemake \
    pandas \
    trimmomatic \
    bowtie2 \
    samtools \
    htseq \
    && micromamba clean --all --yes

WORKDIR /work
COPY --chown=$MAMBA_USER:$MAMBA_USER Go_bacteriaRNake_paired_V4.smk /opt/rnake/Go_bacteriaRNake_paired_V4.smk

ENTRYPOINT ["snakemake", "-s", "/opt/rnake/Go_bacteriaRNake_paired_V4.smk"]
CMD ["--help"]
