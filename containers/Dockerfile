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
    && Rscript -e 'expected <- c(xfun="0.60", lattice="0.23.1", yaml="2.3.12", highr="0.12", evaluate="1.0.5", Rcpp="1.1.2", locfit="1.5.9.12", limma="3.54.2", knitr="1.51", edgeR="3.40.2"); actual <- vapply(names(expected), function(p) as.character(packageVersion(p)), character(1)); stopifnot(identical(actual, expected)); stopifnot(as.character(packageVersion("metaSPARSim")) == "1.1.2"); d <- packageDescription("metaSPARSim"); stopifnot(d$RemoteSha == "37b5931ab1ed6836ded1b605ab23e41824d0f91e")'

# ------------------------------------------------------------
# Install this exact MHASS source tree
# ------------------------------------------------------------

WORKDIR /opt/mhass

COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md ./
COPY --chown=$MAMBA_USER:$MAMBA_USER mhass ./mhass
COPY --chown=$MAMBA_USER:$MAMBA_USER tests ./tests

RUN python -m pip install \
        --no-deps \
        --no-build-isolation \
        .

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
