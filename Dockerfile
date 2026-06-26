FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

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

# Install Python wrapper (python-wrapper branch is a separate deliverable)
RUN git clone -b python-wrapper https://github.com/KidneyRegeneration/DevKidCC /opt/DevKidCC_python
WORKDIR /opt/DevKidCC_python
RUN pip install --no-deps -e .

# Install scPred from source
RUN git clone https://github.com/powellgenomicslab/scPred /opt/scPred
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

# Copy run scripts to a convenient location
COPY docker/run_dkcc.R /opt/run_dkcc.R
COPY docker/run_dkcc_batch.sh /opt/run_dkcc_batch.sh
RUN chmod +x /opt/run_dkcc.R /opt/run_dkcc_batch.sh

WORKDIR /data
CMD ["/bin/bash"]
