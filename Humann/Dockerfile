FROM mambaorg/micromamba:1.5.10

LABEL author="Heekuk Park"
LABEL project="Humann"

ENV ENV_NAME=humann

RUN micromamba create -y -n ${ENV_NAME} \
    -c conda-forge -c bioconda \
    python=3.12 \
    snakemake=8.25.3 \
    humann=3.9 \
    metaphlan=4.1.0 \
    'bowtie2>=2.5,<2.6' \
    libstdcxx-ng \
    libgcc-ng \
    diamond \
    biom-format \
    pandas \
    scipy \
    scikit-learn \
    pigz \
    coreutils \
    pip \
    && BOWTIE_PKG_BIN="$(dirname "$(find /opt/conda/pkgs -path '*/bowtie2-*/bin/bowtie2-build' | sort | tail -n 1)")" \
    && for exe in bowtie2 bowtie2-build bowtie2-align-s bowtie2-align-l bowtie2-build-s bowtie2-build-l bowtie2-inspect bowtie2-inspect-s bowtie2-inspect-l; do \
         if [ -x "$BOWTIE_PKG_BIN/$exe" ]; then \
           ln -sf "$BOWTIE_PKG_BIN/$exe" "/opt/conda/envs/${ENV_NAME}/bin/$exe"; \
         fi; \
       done \
    && micromamba run -n ${ENV_NAME} pip install --no-cache-dir MUSiCC==1.0.4 \
    && micromamba run -n ${ENV_NAME} python -c "import humann; print(humann.__file__)" \
    && micromamba run -n ${ENV_NAME} humann --version \
    && micromamba run -n ${ENV_NAME} python -c "import metaphlan; print(metaphlan.__file__)" \
    && ls -l "/opt/conda/envs/${ENV_NAME}/bin/bowtie2-build" \
    && LD_LIBRARY_PATH="/opt/conda/envs/${ENV_NAME}/lib:$LD_LIBRARY_PATH" PATH="/opt/conda/envs/${ENV_NAME}/bin:$PATH" /opt/conda/envs/${ENV_NAME}/bin/bowtie2-build --version | head -n 1 \
    && micromamba clean --all --yes

# MetaPhlAn 4.x exits non-zero for --version in some environments,
# which breaks HuMAnN 3.x's subprocess.check_output version check.
# Replace bash wrapper with a Python wrapper (avoids /bin/bash dependency).
RUN REAL=/opt/conda/envs/${ENV_NAME}/bin/metaphlan.real && \
    LINK=/opt/conda/envs/${ENV_NAME}/bin/metaphlan && \
    PY=/opt/conda/envs/${ENV_NAME}/bin/python3 && \
    ([ -f "$REAL" ] || mv "$LINK" "$REAL") && \
    printf '#!%s\nimport sys,os,subprocess\n_real=os.path.join(os.path.dirname(os.path.abspath(__file__)),"metaphlan.real")\nif len(sys.argv)>1 and sys.argv[1]=="--version":\n r=subprocess.run([_real]+sys.argv[1:],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)\n sys.stdout.buffer.write(r.stdout)\n sys.exit(0)\nos.execv(_real,[_real]+sys.argv[1:])\n' \
        "$PY" > "$LINK" && \
    chmod +x "$LINK"

ENV MAMBA_DEFAULT_ENV=${ENV_NAME}
ENV PATH=/opt/conda/envs/${ENV_NAME}/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/conda/envs/${ENV_NAME}/lib:$LD_LIBRARY_PATH

WORKDIR /work
CMD ["bash"]
