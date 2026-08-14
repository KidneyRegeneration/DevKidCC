# DevKidCC — Python Wrapper

**Dev**eloping **Kid**ney **C**ell **C**lassifier, callable from Python.

Classify cells in an `AnnData` object with the R DevKidCC package, without
leaving a Scanpy workflow.

[![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![R Version](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This branch (`python-wrapper`) of the DevKidCC repository holds the Python
package. The R package itself lives on `main`; the container that ships both is
built from `containerise`.

## How it works

The wrapper hands the counts matrix to R over a **subprocess**, not rpy2:

```
AnnData --> counts.csv + obs.csv --> Rscript run_dkcc.R --> results.csv --> .obs
```

rpy2 and reticulate segfault when both are loaded into one process, and Seurat
pulls reticulate in. Talking to `Rscript` over a pipe sidesteps that entirely —
at the cost of a CSV round trip, which is the reason for the gene projection
described below.

Classification is hierarchical, and it is `DKCC()` in the R package — not this
wrapper — that implements it:

1. **Stage 1** assigns a broad lineage (Nephron, Stroma, UrEp, NPC, Endo, …)
2. **Stage 2/3** refine within that lineage
3. Cells left `unassigned` are optionally rescued by a KNN vote over the UMAP
   embedding (`knn_smooth`)

## Requirements

**Python** ≥ 3.8, plus `numpy`, `pandas`, `scipy`, `anndata`, `scanpy`
(`pip install -r requirements.txt`). No rpy2, no anndata2ri.

**R** ≥ 4.0 with `Seurat` (v5), `scPred`, and `DevKidCC` **≥ 0.5.1** — earlier
releases have no `knn.iter` parameter and `run_dkcc.R` will refuse to run
against them. Nothing is auto-installed; the classifier checks at construction
and tells you what is missing.

```r
install.packages("Seurat")
remotes::install_github("powellgenomicslab/scPred")
remotes::install_github("KidneyRegeneration/DevKidCC")
```

## Installation

```bash
git clone -b python-wrapper https://github.com/KidneyRegeneration/DevKidCC
cd DevKidCC
pip install .
```

Or skip the R setup entirely and use the container, which carries R, Seurat,
scPred, DevKidCC and this wrapper:

```bash
singularity pull dkcc.sif docker://ghcr.io/kidneyregeneration/dkcc:latest
```

## Quick start

```python
import scanpy as sc
from devkidcc import classify_kidney_cells

adata = sc.read_h5ad("kidney_organoid.h5ad")
adata = classify_kidney_cells(adata)

print(adata.obs[["LineageID", "DKCC"]].value_counts())
```

Input should be **raw counts** with HGNC gene symbols. Strip any genome prefix
first — Cellranger multi-genome output writes `GRCh38_GAPDH`, which matches
nothing in the reference.

## Usage

### Turning KNN smoothing off

The unassigned-cell rescue is on by default. For benchmarking against raw scPred
assignments:

```python
adata_raw = classify_kidney_cells(adata, knn_smooth=False)
```

`knn_iter` takes an explicit iteration count if you want one (it maps straight
onto `DKCC()`'s `knn.iter`; `0` disables smoothing, and is what `knn_smooth=False`
sets).

### Reusing the classifier

Construction runs the R dependency checks, so build it once for a batch:

```python
from devkidcc import DevKidCCClassifier

classifier = DevKidCCClassifier(verbose=False)
for path in samples:
    adata = sc.read_h5ad(path)
    adata = classifier.classify(adata)
    adata.write_h5ad(path.replace(".h5ad", "_DKCC.h5ad"))
```

### Large datasets

The CSV handoff is the memory bottleneck: R loads the whole file before it
builds a Seurat object. Two things keep that in hand.

**The gene projection.** Only DevKidCC reference genes are written — everything
else is discarded by scPred anyway. A 56k × 40k matrix goes from a ~12 GB CSV
and a ~25 GB Seurat object to roughly ~3 GB and ~6 GB. The returned AnnData
still carries every gene you passed in; the projection applies only to what
crosses the process boundary.

Per-cell library sizes are computed over the **full** matrix and sent across
with the metadata, because Seurat's `LogNormalize` otherwise divides by the sum
over the projected genes alone — which moves every normalised value, and every
scPred score with it. On Howden 2019 that was the difference between LineageID
κ 0.66 and κ 0.9989 against unprojected labels. So the projection is a memory
optimisation and nothing more, which is the only thing it should be.

The list ships as `devkidcc/data/reference_genes.txt` (the union of feature
loadings across all seven scPred models). Override it with `$DEVKIDCC_REF_GENES`,
or regenerate it with `Rscript scripts/export_reference_genes.R`.

**Chunking.** Above ~50k cells, classify in chunks and concatenate:

```python
import anndata as ad

chunks = [classifier.classify(adata[i:i + 10_000].copy())
          for i in range(0, adata.n_obs, 10_000)]
classified = ad.concat(chunks)
```

Also point `TMPDIR` somewhere with room — the counts CSV lands there, and a
tmpfs `/tmp` will run out.

## Output

| Column | Description | Example values |
|--------|-------------|----------------|
| `LineageID` | Stage 1 lineage | `Nephron`, `Stroma`, `UrEp`, `NPC`, `Endo`, `unassigned` |
| `DKCC` | Refined cell type | `Podocyte`, `PT`, `EDT`, `SPC`, `NPC-like` |
| `LineageID_max` | Stage 1 confidence | `0.0`–`1.0` |

Per-model scPred probability columns are preserved alongside these.

## Troubleshooting

**"Installed DevKidCC::DKCC() has no knn.iter parameter"** — the R package is
older than 0.5.1. Reinstall from GitHub.

**"None of the input genes are DevKidCC reference genes"** — `var_names` are not
HGNC symbols, or carry a genome prefix. Check `adata.var_names[:10]`.

**"DevKidCC reference gene list not found"** — the packaged list is missing (an
incomplete install). Regenerate it with `scripts/export_reference_genes.R`, or
point `$DEVKIDCC_REF_GENES` at a copy. The wrapper deliberately raises here
rather than falling back to the full matrix, which is what runs a machine out of
memory.

**Killed with no error** — the OOM killer. Reduce the chunk size and check
`TMPDIR`.

## Development

```bash
pip install pytest
pytest                      # fast tests
pytest -m slow -s           # end-to-end; needs R and real data
```

The slow tests are opt-in because they call R. `tests/test_regression_howden.py`
compares a fresh run against stored labels and takes its input from
`$DEVKIDCC_REGRESSION_H5AD` — see its docstring.

## Citation

**Wilson et al., 2022, Genome Medicine** — "Integrated single-cell genomics
reveals the landscape of epithelial, stromal and vascular development in human
fetal kidney"
https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-022-01023-z

## Related

- R DevKidCC package: https://github.com/KidneyRegeneration/DevKidCC
- Container image: https://ghcr.io/kidneyregeneration/dkcc
- Seurat: https://satijalab.org/seurat/
- Scanpy: https://scanpy.readthedocs.io/

## Contact

- **Author**: Sean Wilson — sean.wilson@sund.ku.dk
- **Issues**: https://github.com/KidneyRegeneration/DevKidCC/issues
