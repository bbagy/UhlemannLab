FROM mambaorg/micromamba:1.5.8

# Workdir
WORKDIR /work

# Install tools (conda-forge + bioconda)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    cutadapt \
    minimap2 \
    seqkit \
    vsearch \
    samtools \
    racon \
    blast \
    nanofilt \
    nanoplot \
    python=3.11 \
    openpyxl \
    && micromamba clean --all --yes

# Default shell
SHELL ["/bin/bash", "-lc"]

# Copy pipeline scripts (optional)
# COPY scripts/ /usr/local/bin/

CMD ["bash"]
