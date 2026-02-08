# DevKidCC - Python Wrapper

**Dev**eloping **Kid**ney **C**ell **C**lassifier - Python interface using rpy2

A Python wrapper for the R DevKidCC package, enabling seamless kidney cell classification in Python/Scanpy workflows.

[![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![R Version](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This package provides a Python interface to the R DevKidCC package, allowing you to:
- ✅ Use DevKidCC directly from Python/Scanpy workflows
- ✅ Automatically convert between AnnData and Seurat formats
- ✅ Auto-install R dependencies
- ✅ Get the same results as the R package

**Classification approach:**
1. **Tier 1**: Broad lineage (Nephron, Stroma, Immune, etc.)
2. **Tier 2**: Detailed cell types within each lineage

## Requirements

### Python Requirements
- Python >= 3.8
- scanpy >= 1.9.0
- anndata >= 0.8.0
- rpy2 >= 3.5.0
- anndata2ri >= 1.1.0

### R Requirements
- R >= 4.0
- The following R packages (auto-installed if missing):
  - Seurat
  - SeuratDisk
  - DevKidCC

## Installation

### Step 1: Install R

If you don't have R installed:

**macOS:**
```bash
brew install r
```

**Ubuntu/Debian:**
```bash
sudo apt-get install r-base r-base-dev
```

**Windows:**
Download from [CRAN](https://cran.r-project.org/bin/windows/base/)

### Step 2: Install Python Package

```bash
pip install devkidcc
```

Or from source:
```bash
git clone https://github.com/KidneyRegeneration/DevKidCC-python
cd DevKidCC-python
pip install -e .
```

### Step 3: Verify Installation

```python
python -c "import devkidcc; print('DevKidCC wrapper installed successfully!')"
```

## Quick Start

```python
import scanpy as sc
from devkidcc import classify_kidney_cells

# Load your kidney organoid or tissue data
adata = sc.read_h5ad("kidney_organoid.h5ad")

# Classify in one line! 
# (R packages will be auto-installed on first run)
adata = classify_kidney_cells(adata)

# View results
print(adata.obs[['LineageID', 'DKCC']].value_counts())
```

That's it! The first run will automatically:
1. Check for required R packages
2. Install them if missing
3. Convert your data to Seurat format
4. Run DevKidCC classification
5. Convert results back to AnnData

## Usage

### Basic Classification

```python
from devkidcc import DevKidCCClassifier

# Create classifier (checks/installs R packages)
classifier = DevKidCCClassifier(verbose=True)

# Classify your data
adata = classifier.classify(adata)

# Results are in adata.obs:
# - 'LineageID': Broad cell lineage
# - 'DKCC': Detailed cell type
```

### Recommended Workflow with QC

```python
import scanpy as sc
from devkidcc import classify_kidney_cells

# Load data
adata = sc.read_h5ad("kidney_data.h5ad")

# Quality control (recommended)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Filter poor quality cells
adata = adata[adata.obs.pct_counts_mt < 20, :]

# Classify
adata = classify_kidney_cells(adata)

# Save
adata.write_h5ad("kidney_classified.h5ad")
```

### Batch Processing

```python
from devkidcc import DevKidCCClassifier

# Create classifier once (R packages loaded once)
classifier = DevKidCCClassifier(verbose=False)

samples = ['sample1.h5ad', 'sample2.h5ad', 'sample3.h5ad']
results = []

for sample_file in samples:
    adata = sc.read_h5ad(sample_file)
    adata = classifier.classify(adata)
    adata.obs['sample'] = sample_file
    results.append(adata)

# Combine
combined = sc.concat(results)
```

### Visualization

```python
import scanpy as sc

# Compute UMAP if needed
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Plot lineages
sc.pl.umap(adata, color='LineageID', legend_loc='right margin')

# Plot detailed cell types
sc.pl.umap(adata, color='DKCC', legend_loc='on data', 
           legend_fontsize=6)

# Cell type proportions
import matplotlib.pyplot as plt
adata.obs['DKCC'].value_counts().plot(kind='barh')
plt.xlabel('Number of Cells')
plt.tight_layout()
plt.show()
```

## How It Works

1. **AnnData → Seurat**: Your AnnData object is saved as h5ad, then converted to Seurat's h5seurat format
2. **R Classification**: The R DevKidCC package classifies cells in Seurat
3. **Seurat → AnnData**: Results are converted back to h5ad and loaded into Python
4. **Cleanup**: Temporary files are automatically removed

```
Python (AnnData) → h5ad → h5seurat → R DevKidCC → h5seurat → h5ad → Python (AnnData)
```

## Performance

- **First run**: 2-5 minutes (R package installation)
- **Subsequent runs**: 1-3 minutes for 10,000 cells
- **Conversion overhead**: ~30 seconds for 10,000 cells

## Comparison with Pure Python Approach

### Advantages
- ✅ Uses the exact same R models (guaranteed identical results)
- ✅ No model retraining needed
- ✅ Automatic R package management
- ✅ Quick to set up (<5 minutes)
- ✅ Always stays in sync with R package updates

### Disadvantages
- ❌ Requires R installation
- ❌ Slower due to format conversions
- ❌ More complex dependencies

## Troubleshooting

### "R is not installed"

Install R from [CRAN](https://cran.r-project.org/) or use your package manager:
```bash
# macOS
brew install r

# Ubuntu
sudo apt-get install r-base
```

### "Unable to install R package"

Install manually in R:
```r
# Open R console
install.packages("Seurat")
install.packages("devtools")
devtools::install_github("mojaveazure/seurat-disk")
devtools::install_github("KidneyRegeneration/DevKidCC")
```

Then try again in Python:
```python
classifier = DevKidCCClassifier(install_deps=False)
```

### "Gene names don't match"

Ensure your gene names are in the correct format (usually HGNC symbols):
```python
# Check current format
print(adata.var_names[:10])

# Convert if needed
adata.var_names = adata.var_names.str.upper()
```

### "Conversion failed"

Check SeuratDisk installation:
```r
# In R
library(Seurat)
library(SeuratDisk)

# If error, reinstall:
devtools::install_github("mojaveazure/seurat-disk", force=TRUE)
```

### Memory Issues

For very large datasets (>100k cells), consider:
```python
# Process in chunks
from sklearn.model_selection import KFold

kf = KFold(n_splits=5)
results = []

for train_idx, test_idx in kf.split(adata.obs):
    subset = adata[test_idx].copy()
    subset = classifier.classify(subset)
    results.append(subset)

adata_classified = sc.concat(results)
```

### rpy2 Installation Issues

**macOS:** If you get compiler errors:
```bash
brew install pkg-config
export PKG_CONFIG_PATH="/usr/local/lib/pkgconfig"
pip install rpy2
```

**Linux:** Install R development headers:
```bash
sudo apt-get install r-base-dev
pip install rpy2
```

## Output Format

Classification adds these columns to `adata.obs`:

| Column | Description | Example Values |
|--------|-------------|----------------|
| `LineageID` | Broad lineage | "Nephron", "Stroma", "Immune" |
| `DKCC` | Detailed cell type | "Podocyte", "Proximal Tubule" |

Additional metadata columns from the R package are also preserved.

## Citation

If you use DevKidCC in your research, please cite:

**Wilson et al., 2022, Genome Medicine**  
"Integrated single-cell genomics reveals the landscape of epithelial, stromal and vascular development in human fetal kidney"  
https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-022-01023-z

## Examples

See the [examples](examples/) directory for:
- `example_usage.py` - Comprehensive examples
- `tutorial.ipynb` - Jupyter notebook tutorial
- `batch_processing.py` - Process multiple samples

## Development

### Running Tests

```bash
pip install pytest
pytest tests/
```

### Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests
4. Submit a pull request

## Related Projects

- **R DevKidCC Package**: https://github.com/KidneyRegeneration/DevKidCC
- **Seurat**: https://satijalab.org/seurat/
- **Scanpy**: https://scanpy.readthedocs.io/

## License

MIT License - see [LICENSE](LICENSE) file

## Contact

- **Author**: Sean Wilson
- **Email**: sean.wilson@mcri.edu.au
- **Issues**: https://github.com/KidneyRegeneration/DevKidCC-python/issues
- **Original R Package**: https://github.com/KidneyRegeneration/DevKidCC

## Acknowledgments

- Original DevKidCC R package by Sean Wilson
- Built using rpy2, Scanpy, and Seurat
- Thanks to the kidney research community