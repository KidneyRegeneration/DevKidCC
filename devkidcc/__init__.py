"""
DevKidCC: Developing Kidney Cell Classifier (Python Wrapper)
=============================================================

A Python wrapper for the R DevKidCC package, enabling kidney cell
classification in Python using AnnData/Scanpy workflows.

This package uses subprocess to call R scripts, avoiding rpy2/reticulate
conflicts. It automatically handles conversions between Python (AnnData)
and R (Seurat) data formats.

Basic Usage
-----------
>>> import scanpy as sc
>>> from devkidcc import classify_kidney_cells
>>>
>>> # Load your data
>>> adata = sc.read_h5ad("kidney_data.h5ad")
>>>
>>> # Classify cells (automatic R package setup)
>>> adata = classify_kidney_cells(adata)
>>>
>>> # View results
>>> print(adata.obs[['LineageID', 'DKCC']].value_counts())

Advanced Usage
--------------
>>> from devkidcc import DevKidCCClassifier
>>>
>>> # Create classifier (installs R dependencies if needed)
>>> classifier = DevKidCCClassifier(verbose=True)
>>>
>>> # Classify
>>> adata = classifier.classify(adata)

Requirements
------------
- Python >= 3.8
- R >= 4.0
- R packages: Seurat, SeuratDisk, DevKidCC (auto-installed)

References
----------
Wilson et al., 2022, Genome Medicine
https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-022-01023-z
"""

__version__ = "0.5.0"
__author__ = "Sean Wilson"
__email__ = "sean.wilson@mcri.edu.au"

from .classifier_subprocess import DevKidCCClassifier, classify_kidney_cells
from .dkcc_v2_plotting_functions import (
    plot_scpred_scores_distribution,
    plot_scpred_umap_py,
    create_custom_cmap
)

__all__ = [
    "DevKidCCClassifier",
    "classify_kidney_cells",
    "plot_scpred_scores_distribution",
    "plot_scpred_umap_py",
    "create_custom_cmap",
]