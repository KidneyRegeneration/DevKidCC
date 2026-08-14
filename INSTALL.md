# Installation

Three ways in, in descending order of how much of your afternoon they cost.

| Route | Use when |
|---|---|
| [Container](#container-recommended) | You want it working now, or you are on an HPC |
| [Existing R install](#existing-r-install) | You already run Seurat locally |
| [From scratch](#from-scratch) | Neither of the above |

The Python side is small — `numpy`, `pandas`, `scipy`, `anndata`, `scanpy`. All
the difficulty is in the R side, and the container exists to remove it.

---

## Container (recommended)

The published image carries R, Seurat, scPred, DevKidCC and this wrapper:

```bash
# Singularity / Apptainer (HPC)
singularity pull dkcc.sif docker://ghcr.io/kidneyregeneration/dkcc:latest
singularity exec dkcc.sif Rscript /opt/run_dkcc.R --input sample.h5ad --output sample_DKCC.h5ad

# Docker
docker run -v "$PWD:/data" ghcr.io/kidneyregeneration/dkcc:latest \
    Rscript /opt/run_dkcc.R --input /data/sample.h5ad --output /data/sample_DKCC.h5ad
```

`run_dkcc.sh` in the `containerise` branch wraps both, including SLURM
submission. Pull the SIF on a login node — compute nodes typically have no
outbound network.

---

## Existing R install

```bash
git clone -b python-wrapper https://github.com/KidneyRegeneration/DevKidCC
cd DevKidCC
pip install .
```

Then in R:

```r
install.packages(c("Seurat", "remotes"))
remotes::install_github("powellgenomicslab/scPred")
remotes::install_github("KidneyRegeneration/DevKidCC")
```

Check what you ended up with — the wrapper requires DevKidCC ≥ 0.5.1 and will
refuse to run against anything older:

```r
packageVersion("DevKidCC")
"knn.iter" %in% names(formals(DevKidCC::DKCC))   # must be TRUE
```

---

## From scratch

### 1. R

```bash
# Ubuntu / Debian
sudo apt-get install r-base r-base-dev

# Fedora / RHEL
sudo dnf install R R-devel

# macOS
brew install r
```

Windows: download from [CRAN](https://cran.r-project.org/bin/windows/base/).

Seurat needs system libraries for its compiled dependencies. On Ubuntu:

```bash
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
```

`environment.yml` on the `containerise` branch is the authoritative list — it is
what the image is built from.

### 2. Python

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
pip install .
```

### 3. R packages

As under [Existing R install](#existing-r-install). Seurat takes a while to
compile the first time.

---

## Verify

```bash
python -c "import devkidcc; print(devkidcc.__version__)"
python -c "from devkidcc import DevKidCCClassifier; DevKidCCClassifier()"
```

The second line runs the R dependency checks and names anything missing. Then
the fast test suite, which needs no R:

```bash
pip install pytest
pytest
```

And end-to-end, which does:

```bash
pytest -m slow -s
```

---

## Troubleshooting

**`Rscript not found`** — R is not on `PATH`. Pass the path explicitly:

```python
DevKidCCClassifier(rscript_path="/usr/local/bin/Rscript")
```

**`Required R package 'DevKidCC' not found`** — installed for a different R
version than the one `Rscript` resolves to. Compare `R --version` with
`.libPaths()` inside `R`.

**`Installed DevKidCC::DKCC() has no knn.iter parameter`** — the R package
predates 0.5.1. Reinstall from GitHub.

**Seurat fails to compile** — missing system libraries; see step 1. This is the
single most common reason a from-scratch install fails, and the reason the
container exists.

**Killed, with no error message** — the OOM killer, during the CSV handoff.
Classify in chunks and put `TMPDIR` somewhere with room; `/tmp` is a
memory-backed tmpfs on many systems, so a large counts CSV written there
consumes the RAM you are trying to conserve.

**`No such file or directory: run_dkcc.R`** — an incomplete install; the R
script is package data. Reinstall rather than copying files around.

---

## Uninstall

```bash
pip uninstall devkidcc
R -e "remove.packages(c('DevKidCC', 'scPred'))"
```

---

## Getting help

Open an issue at
https://github.com/KidneyRegeneration/DevKidCC/issues with your OS, `R --version`,
`python --version`, `packageVersion("DevKidCC")`, and the full error.
