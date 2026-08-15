
# DevKidCC

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![Seurat v5](https://img.shields.io/badge/Seurat-v5-blue.svg)](https://satijalab.org/seurat/)
[![Version](https://img.shields.io/badge/version-0.5.1-blue.svg)](https://github.com/KidneyRegeneration/DevKidCC)
<!-- badges: end -->

**DevKidCC** (*Dev*eloping *Kid*ney *C*ell *C*lassifier) is a tool that will classify single cell kidney data, both human tissue and human stem cell derived organoids. There is no pre-processing required, although we do recommend filtering out poor quality cells for most accurate representation of cell proportions.

**Version 0.5.0** introduces a simplified NPC (Nephron Progenitor Cell) refinement approach for improved Seurat v5 compatibility. NPC cells are now classified as "NPC-like" without PAX2-based clustering refinement. This change resolves data frame compatibility issues while maintaining classification accuracy for all other cell types. See [NEWS.md](NEWS.md) for full details.

**Backward Compatibility Note:** If you have existing analyses using v0.4.0, NPC cell classifications may differ in v0.5.0. Cells previously labeled as "NPC" based on PAX2 expression clustering are now uniformly labeled as "NPC-like". All other cell type classifications remain unchanged. This change improves Seurat v5 compatibility and eliminates data frame errors but removes NPC subtyping granularity.

bioRxiv paper: [Wilson et al. 2021](https://doi.org/10.1101/2021.01.20.427346) <br>
Genome Medicine: [Wilson et al., 2022](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-022-01023-z)

<br>
<img src="./Model_Graphic.png">
<br>


## Installation

### Option 1 — Container (no R or Python installation required)

A pre-built image is published to the GitHub Container Registry on every commit to the release branches. It carries R, Seurat, scPred and DevKidCC alongside Python, anndata and the `devkidcc` wrapper, so **either language can drive the classifier**.

```bash
docker pull ghcr.io/kidneyregeneration/dkcc:latest
```

Check the image is sound before trusting it with your data. This pushes a synthetic matrix through both entry points and takes about two minutes; it needs no input files and no network:

```bash
docker run --rm ghcr.io/kidneyregeneration/dkcc:latest python /opt/smoke_test.py
# ... PASS: both the R and Python entry points classified the input.
```

**R entry point** — mount your data directory to `/data` inside the container:

```bash
# Run on a single file
docker run --rm -v /path/to/your/data:/data \
    ghcr.io/kidneyregeneration/dkcc:latest \
    Rscript /opt/run_dkcc.R \
        --input  /data/sample.h5ad \
        --output /data/sample_DKCC.h5ad \
        --format h5ad

# Batch process all supported files in a directory
docker run --rm -v /path/to/your/data:/data \
    ghcr.io/kidneyregeneration/dkcc:latest \
    bash /opt/run_dkcc_batch.sh -i /data -o h5ad
```

Supported input formats: `.h5ad`, `.h5`, `.h5seurat`, `.rds`, `.RData`  
Supported output formats: `h5ad`, `rds`

**Python entry point** — for calling DevKidCC from a scanpy workflow. The wrapper shells out to R for you; `DKCC` and `LineageID` come back as `.obs` columns on your AnnData, at its original gene width:

```bash
docker run --rm -v /path/to/your/data:/data \
    ghcr.io/kidneyregeneration/dkcc:latest python -c "
import anndata as ad, devkidcc
adata = ad.read_h5ad('/data/sample.h5ad')
adata = devkidcc.classify_kidney_cells(adata)
print(adata.obs['LineageID'].value_counts())
adata.write_h5ad('/data/sample_DKCC.h5ad')
"
```

Pass raw counts, not log-normalised values — the models were fit against counts, and normalisation happens inside.

**Singularity / Apptainer**, for HPC systems where Docker is unavailable:

```bash
singularity pull dkcc.sif docker://ghcr.io/kidneyregeneration/dkcc:latest
singularity exec --bind /path/to/data:/data dkcc.sif \
    Rscript /opt/run_dkcc.R --input /data/sample.h5ad --output /data/sample_DKCC.h5ad
```

`run_dkcc.sh` in this repository wraps that call for SLURM clusters — see `HPC_TESTING.md`. Note that Singularity bind-mounts your home directory by default; the image sets `PYTHONNOUSERSITE=1` so a `~/.local` Python install on the host cannot shadow the versions inside it.

The image is built for `linux/amd64`. On Apple Silicon it runs under emulation, slowly.

### Option 2 — R package (devtools)

You can install **DevKidCC** from this repository using devtools: 

``` r
# prerequite scPred package
devtools::install_github("powellgenomicslab/scPred")
# DevKidCC itself
devtools::install_github("KidneyRegeneration/DevKidCC", ref = "main")
```
Expected installation time is under 5 minutes for each package.

This package has been successfully tested on both Windows and Linux systems.

## Standard Workflow

Check out the full vignette which includes details on using the visualisation functions [here](https://kidneyregeneration.github.io/DevKidCC/index.html)

### Running DevKidCC

DevKidCC operates on single cell data as a Seurat object. The simplest workflow to use DevKidCC is:

``` r
library(DevKidCC)
# read in seurat object
organoid <- DKCC(organoid)  # use dataset 'organoid' included in package
```
This will cause a number of additional metadata to be added to the object. The first tier is labelled "LineageID" while the complete annotation is under "DKCC"


### Loading organoid gene expression database

The database is stored as an rda file and can be downloaded at the following link:
https://drive.google.com/file/d/1wh551HvecgszizE8FCsXRD5CiQVYm3K6/view?usp=sharing

Save this file in a data folder at the top level of your project. This will allow it to be loaded by the DotPlotCompare function when `compare.to.organoids = T`



Any issues submit a pull request or contact me at sean.wilson@mcri.edu.au

