FROM mambaorg/micromamba:1.5.10

LABEL author="Heekuk Park"
LABEL project="RNake"

ENV ENV_NAME=rnake

RUN micromamba create -y -n ${ENV_NAME} \
    -c conda-forge -c bioconda \
    python=3.11 \                        
    snakemake \                          
    pandas \                             
    trimmomatic \                        
    bowtie2 \                            
    samtools \                           
    htseq \                              
    && micromamba clean --all --yes   

ENV PATH=/opt/conda/envs/${ENV_NAME}/bin:$PATH
ENV CONDA_DEFAULT_ENV=${ENV_NAME}

WORKDIR /work
CMD ["bash"]