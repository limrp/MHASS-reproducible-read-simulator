# Pinned linux/amd64 manifest for mambaorg/micromamba:2.9.0
FROM mambaorg/micromamba:2.9.0@sha256:5681ae3caa12844c41d31e1d70f636fbaba31b01bdfd1b5b521742d7e490614d

# ------------------------------------------------------------
# Recreate the exact known-good Conda portion of the environment
# ------------------------------------------------------------

COPY --chown=$MAMBA_USER:$MAMBA_USER \
    containers/environment-linux-64.lock \
    /tmp/environment-linux-64.lock

RUN micromamba install \
        --yes \
        --name base \
        --file /tmp/environment-linux-64.lock \
    && micromamba clean --all --yes

# Install Samtools in an isolated Conda environment.
COPY --chown=$MAMBA_USER:$MAMBA_USER \
    containers/samtools-linux-64.lock \
    /tmp/samtools-linux-64.lock

RUN micromamba create \
        --yes \
        --name samtools \
        --file /tmp/samtools-linux-64.lock \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/envs/samtools/bin:${PATH}"


# ------------------------------------------------------------
# Make the base Conda environment available to Docker RUN steps
# ------------------------------------------------------------

ARG MAMBA_DOCKERFILE_ACTIVATE=1


# ------------------------------------------------------------
# Recreate required pip-installed packages
# ------------------------------------------------------------

COPY --chown=$MAMBA_USER:$MAMBA_USER \
    containers/requirements-pip.txt \
    /tmp/requirements-pip.txt

RUN python -m pip install \
        --no-cache-dir \
        --no-deps \
        -r /tmp/requirements-pip.txt

# ------------------------------------------------------------
# Install the exact BiocManager version required by get_counts.R
# ------------------------------------------------------------

RUN curl \
        --fail \
        --location \
        --show-error \
        https://cran.r-project.org/src/contrib/BiocManager_1.30.27.tar.gz \
        --output /tmp/BiocManager_1.30.27.tar.gz \
    && echo "acbfdcf09602c8279085556ca54531c2ada0ac3c1775d54ce2c3b9a3dc03fbb5  /tmp/BiocManager_1.30.27.tar.gz" \
        | sha256sum --check \
    && R CMD INSTALL /tmp/BiocManager_1.30.27.tar.gz \
    && Rscript -e 'stopifnot(as.character(packageVersion("BiocManager")) == "1.30.27")' \
    && rm /tmp/BiocManager_1.30.27.tar.gz

# ------------------------------------------------------------
# Install the audited metaSPARSim revision
# ------------------------------------------------------------

RUN Rscript -e 'remotes::install_gitlab("sysbiobig/metaSPARSim", ref="37b5931ab1ed6836ded1b605ab23e41824d0f91e", upgrade="never")' \
    && Rscript -e 'stopifnot(as.character(packageVersion("metaSPARSim")) == "1.1.2")'


# ------------------------------------------------------------
# Install this exact MHASS source tree
# ------------------------------------------------------------

WORKDIR /opt/mhass

COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md ./
COPY --chown=$MAMBA_USER:$MAMBA_USER mhass ./mhass
COPY --chown=$MAMBA_USER:$MAMBA_USER tests ./tests

RUN python -m pip install --no-deps .


# ------------------------------------------------------------
# PBSIM resource location
# ------------------------------------------------------------

ENV PBSIM_QSHMM=/opt/conda/data/QSHMM-RSII.model


# ------------------------------------------------------------
# Build-time validation
# ------------------------------------------------------------

RUN command -v python \
    && command -v R \
    && command -v pbsim \
    && command -v ccs \
    && command -v samtools \
    && samtools --version \
    && test -f "$PBSIM_QSHMM" \
    && pytest -q tests/

CMD ["mhass", "-h"]
