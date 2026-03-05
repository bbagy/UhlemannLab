  GNU nano 7.2                                                                                         Dockerfile                                                                                                   
FROM mambaorg/micromamba:1.5.10

LABEL author="Heekuk Park"
LABEL project="PF-Snake"

# Environment name
ENV ENV_NAME=pf-snake

# Install bioinformatics tools
RUN micromamba create -y -n ${ENV_NAME} \
    -c conda-forge -c bioconda \
    python=3.10 \
    snakemake \
    bwa \
    samtools \
    bcftools \
    freebayes \
    picard \
    bedtools \
    seqtk \
    fastp \
    gatk4 \
    mafft \
    muscle \
    vcftools \
    tabix \
    htslib \
    pandas \
    numpy \
    multiqc \
    openpyxl \
    && micromamba clean --all --yes

# Activate env in the default shell
ENV PATH=/opt/conda/envs/${ENV_NAME}/bin:$PATH

WORKDIR /work
