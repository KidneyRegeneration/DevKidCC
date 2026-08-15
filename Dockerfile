FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# R emits seven "Setting LC_* failed" warnings per invocation on a locale-less
# image, which buries real output in every log the user sends us.
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# Singularity bind-mounts $HOME by default, which puts the *host's*
# ~/.local/lib/pythonX.Y/site-packages on sys.path ahead of this image's
# environment. An external user with their own numpy or anndata there would
# silently run it instead of the pinned one -- the classic "works for me" HPC
# failure. Ignore user site-packages outright; nothing here installs into it.
ENV PYTHONNOUSERSITE=1

RUN apt-get update && apt-get install -y \
    wget \
    git \
    curl \
    build-essential \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install micromamba
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xvj bin/micromamba \
    && mv bin/micromamba /usr/local/bin/

ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV MAMBA_EXE=/usr/local/bin/micromamba

COPY environment.yml /tmp/environment.yml
RUN micromamba create -y -n devkid -f /tmp/environment.yml

ENV PATH=/opt/micromamba/envs/devkid/bin:$PATH

# Install additional R packages (remotes, scCustomize)
COPY docker/install.R /tmp/install.R
RUN Rscript /tmp/install.R

# Install the Python wrapper, which lives on its own branch of this same repo.
# Non-editable: an editable install leaves the image depending on /opt/DevKidCC_python
# staying present and writable, which buys nothing here.
#
# The layer cache keys on the instruction text, not on what the remote branch
# points at, so building with the branch *name* here would keep serving whatever
# that branch held the first time this image built -- new wrapper commits would
# silently never reach the image. CI therefore resolves the branch to a commit
# SHA and passes it as --build-arg WRAPPER_REF=<sha>, which changes the
# instruction text exactly when the branch moves. The branch name below is only
# a fallback for local builds; pass a SHA or tag to get a reproducible one.
ARG WRAPPER_REF=python-wrapper-v0.5.1
RUN git clone https://github.com/KidneyRegeneration/DevKidCC /opt/DevKidCC_python \
    && git -C /opt/DevKidCC_python checkout --quiet ${WRAPPER_REF}
WORKDIR /opt/DevKidCC_python
RUN pip install --no-deps .

# Install scPred from source. Pinned: scPred's master has no releases, so an
# unpinned clone makes every rebuild of this image a different image.
ARG SCPRED_REF=af5492e778b076e529c20462c92aacd06c75bdc0
RUN git clone https://github.com/powellgenomicslab/scPred /opt/scPred \
    && git -C /opt/scPred checkout --quiet ${SCPRED_REF}
WORKDIR /opt/scPred
RUN Rscript -e 'remotes::install_local("/opt/scPred", upgrade="never")'

# Install DevKidCC R package from the local build context
COPY . /opt/DevKidCC
WORKDIR /opt/DevKidCC
RUN rm -f .Rprofile
RUN Rscript -e 'remotes::install_local(".", upgrade="never")'

# Sanity checks
RUN python --version
RUN python -c "import devkidcc"

# run_dkcc.R reads h5ad through sceasy, which reaches Python's anndata via
# reticulate. Until now this variable was set only by run_dkcc.sh at exec time,
# so reticulate had to go looking whenever the R script was invoked directly --
# which is exactly how the README documents the R entry point, and how anyone
# running `docker run ... Rscript /opt/run_dkcc.R` reaches it. Set it in the
# image so the entry point works on its own terms.
#
# Deliberately placed after the heavy layers: an ENV near the top of the file
# invalidates every layer below it, and this one is only needed at runtime.
ENV RETICULATE_PYTHON=/opt/micromamba/envs/devkid/bin/python

# Copy run scripts to a convenient location
COPY docker/run_dkcc.R /opt/run_dkcc.R
COPY docker/run_dkcc_batch.sh /opt/run_dkcc_batch.sh
COPY docker/smoke_test.py /opt/smoke_test.py
RUN chmod +x /opt/run_dkcc.R /opt/run_dkcc_batch.sh /opt/smoke_test.py

WORKDIR /data
CMD ["/bin/bash"]
