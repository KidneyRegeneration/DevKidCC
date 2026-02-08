"""
DevKidCC Python Wrapper
Provides a Python interface to the R DevKidCC package using rpy2
"""

import warnings
import tempfile
import os
from pathlib import Path
from typing import Optional, Union
import numpy as np
import pandas as pd
import anndata as ad

# RPy2 imports
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter

# Note: pandas2ri and numpy2ri conversions are handled via localconverter
# The activate() method is deprecated in rpy2 >= 3.5
# Conversions are now applied via context managers when needed


class DevKidCCClassifier:
    """
    Python wrapper for the DevKidCC R package.
    
    This class provides a Python interface to classify kidney cells using
    the R DevKidCC package. It handles conversion between AnnData (Python)
    and Seurat (R) objects automatically.
    
    Parameters
    ----------
    verbose : bool, default=True
        Whether to print progress messages
    install_deps : bool, default=True
        Whether to automatically install R dependencies if missing
        
    Examples
    --------
    >>> import scanpy as sc
    >>> from devkidcc import DevKidCCClassifier
    >>> 
    >>> # Load data
    >>> adata = sc.read_h5ad("kidney_organoid.h5ad")
    >>> 
    >>> # Classify
    >>> classifier = DevKidCCClassifier()
    >>> adata = classifier.classify(adata)
    >>> 
    >>> # View results
    >>> print(adata.obs[['LineageID', 'DKCC']].value_counts())
    """
    
    def __init__(self, verbose: bool = True, install_deps: bool = True):
        self.verbose = verbose
        self._r_packages_loaded = False
        
        if self.verbose:
            print("Initializing DevKidCC Python wrapper...")
        
        # Initialize R packages
        self._setup_r_packages(install_deps)
    
    def _setup_r_packages(self, install_deps: bool = True):
        """Setup and load required R packages."""

        # Import R utilities
        utils = importr('utils')
        base = importr('base')

        # Prevent reticulate segfault issue with rpy2
        # Set option to not auto-load reticulate
        ro.r('options(reticulate.autoload = FALSE)')

        # Update reticulate to latest version to avoid segfault with rpy2
        # See: https://github.com/rstudio/reticulate/pull/1188
        if self.verbose:
            print("\n  Updating reticulate to avoid segfault...")
        try:
            utils.install_packages('reticulate', repos='https://cloud.r-project.org')
            if self.verbose:
                print("  [OK] reticulate updated")
        except:
            if self.verbose:
                print("  [WARNING] Could not update reticulate")

        required_packages = {
            'remotes': 'CRAN',
            'Seurat': 'CRAN',
            'scPred': 'github:powellgenomicslab/scPred',
            'SeuratDisk': 'github:mojaveazure/seurat-disk',
            'DevKidCC': 'github:KidneyRegeneration/DevKidCC'
        }
        
        if self.verbose:
            print("\nChecking R package dependencies...")
        
        for pkg, source in required_packages.items():
            try:
                if self.verbose:
                    print(f"  Checking {pkg}...", end=" ")

                # Try to check if package is installed (avoid importr for Seurat to prevent segfault)
                check_installed = ro.r(f'requireNamespace("{pkg}", quietly = TRUE)')[0]

                if not check_installed:
                    raise Exception("Package not installed")

                if self.verbose:
                    print("[OK]")

            except Exception as e:
                if self.verbose:
                    print(f"[X] (not installed)")
                
                if install_deps:
                    if self.verbose:
                        print(f"    Installing {pkg} from {source}...")
                    
                    try:
                        if source == 'CRAN':
                            utils.install_packages(pkg, repos='https://cloud.r-project.org')
                        elif source.startswith('github:'):
                            # Install remotes first if needed
                            try:
                                importr('remotes')
                            except:
                                utils.install_packages('remotes', repos='https://cloud.r-project.org')

                            remotes = importr('remotes')
                            repo = source.replace('github:', '')
                            remotes.install_github(repo)
                        
                        if self.verbose:
                            print(f"    [OK] {pkg} installed successfully")
                    except Exception as install_error:
                        raise RuntimeError(
                            f"Failed to install R package {pkg}: {str(install_error)}\n"
                            f"Please install manually in R: install.packages('{pkg}')"
                        )
                else:
                    raise RuntimeError(
                        f"Required R package {pkg} not found. "
                        f"Please install it in R or set install_deps=True"
                    )
        
        # Load the packages via library() to avoid reticulate segfault with importr()
        ro.r('library(Seurat)')
        ro.r('library(SeuratDisk)')
        ro.r('library(DevKidCC)')

        # Store references (using ro.r to access functions)
        self.r_base = base
        
        self._r_packages_loaded = True
        
        if self.verbose:
            print("\n[OK] All R packages loaded successfully\n")
    
    def _adata_to_seurat(self, adata: ad.AnnData, filename: Optional[str] = None) -> str:
        """
        Convert AnnData to Seurat object via h5ad/h5seurat.
        
        Parameters
        ----------
        adata : anndata.AnnData
            Input AnnData object
        filename : str, optional
            Base filename (without extension). If None, uses temp file.
            
        Returns
        -------
        str
            Path to the h5seurat file
        """
        if filename is None:
            # Create temporary file
            temp_dir = tempfile.mkdtemp()
            filename = os.path.join(temp_dir, "temp_data")
        
        h5ad_path = f"{filename}.h5ad"
        h5seurat_path = f"{filename}.h5seurat"
        
        if self.verbose:
            print(f"Converting AnnData to Seurat...")
            print(f"  Saving h5ad: {h5ad_path}")
        
        # Save AnnData as h5ad
        adata.write_h5ad(h5ad_path)
        
        # Convert to h5seurat using R
        if self.verbose:
            print(f"  Converting to h5seurat...")
        
        ro.r(f'''
            library(Seurat)
            library(SeuratDisk)
            
            # Convert h5ad to h5seurat
            Convert("{h5ad_path}", dest = "h5seurat", overwrite = TRUE)
        ''')
        
        if self.verbose:
            print(f"  [OK] Conversion complete: {h5seurat_path}")
        
        return h5seurat_path
    
    def _seurat_to_adata(self, h5seurat_path: str) -> ad.AnnData:
        """
        Convert Seurat h5seurat file back to AnnData.
        
        Parameters
        ----------
        h5seurat_path : str
            Path to h5seurat file
            
        Returns
        -------
        anndata.AnnData
            Converted AnnData object
        """
        # Convert h5seurat back to h5ad
        base_path = h5seurat_path.replace('.h5seurat', '')
        h5ad_path = f"{base_path}.h5ad"
        
        if self.verbose:
            print(f"\nConverting results back to AnnData...")
        
        ro.r(f'''
            library(Seurat)
            library(SeuratDisk)
            
            # Convert h5seurat to h5ad
            Convert("{h5seurat_path}", dest = "h5ad", overwrite = TRUE)
        ''')
        
        # Read back into Python
        adata = ad.read_h5ad(h5ad_path)
        
        if self.verbose:
            print(f"  [OK] Loaded results back to Python")
        
        return adata
    
    def classify(self, adata: ad.AnnData, copy: bool = True) -> ad.AnnData:
        """
        Classify kidney cells using DevKidCC.
        
        This is the main classification method. It converts the AnnData object
        to Seurat format, runs the R DevKidCC classification, and converts
        the results back to AnnData.
        
        Parameters
        ----------
        adata : anndata.AnnData
            Input single-cell RNA-seq data
        copy : bool, default=True
            Whether to return a copy (recommended)
            
        Returns
        -------
        anndata.AnnData
            Annotated data with classifications in .obs:
            - 'LineageID': Broad lineage category
            - 'DKCC': Detailed cell type annotation
            
        Examples
        --------
        >>> import scanpy as sc
        >>> from devkidcc import DevKidCCClassifier
        >>> 
        >>> adata = sc.read_h5ad("kidney_data.h5ad")
        >>> classifier = DevKidCCClassifier()
        >>> adata = classifier.classify(adata)
        """
        if not self._r_packages_loaded:
            raise RuntimeError("R packages not loaded. Initialization failed.")
        
        if copy:
            adata = adata.copy()
        
        if self.verbose:
            print("="*60)
            print("DevKidCC Classification Pipeline (R backend)")
            print("="*60)
            print(f"Input: {adata.n_obs} cells × {adata.n_vars} genes\n")
        
        # Create temp directory for conversion files
        temp_dir = tempfile.mkdtemp()
        base_path = os.path.join(temp_dir, "data")
        
        try:
            # Convert AnnData to Seurat
            h5seurat_path = self._adata_to_seurat(adata, base_path)
            
            # Run DevKidCC in R
            if self.verbose:
                print("\nRunning DevKidCC classification in R...")
            
            ro.r(f'''
                library(Seurat)
                library(SeuratDisk)
                library(DevKidCC)
                
                # Load Seurat object
                seurat_obj <- LoadH5Seurat("{h5seurat_path}")
                
                # Run DevKidCC classification
                seurat_obj <- DKCC(seurat_obj)
                
                # Save back to h5seurat
                SaveH5Seurat(seurat_obj, "{h5seurat_path}", overwrite = TRUE)
            ''')
            
            if self.verbose:
                print("  [OK] Classification complete")
            
            # Convert back to AnnData
            result = self._seurat_to_adata(h5seurat_path)
            
            # Copy relevant obs columns to original adata
            classification_cols = ['LineageID', 'DKCC']
            for col in classification_cols:
                if col in result.obs.columns:
                    adata.obs[col] = result.obs[col]
            
            # Also copy any score columns if they exist
            score_cols = [c for c in result.obs.columns if 'score' in c.lower() or 'prob' in c.lower()]
            for col in score_cols:
                adata.obs[col] = result.obs[col]
            
            if self.verbose:
                print("\n" + "="*60)
                print("Classification Results")
                print("="*60)
                if 'LineageID' in adata.obs:
                    print(f"\nLineages identified: {adata.obs['LineageID'].nunique()}")
                    print(adata.obs['LineageID'].value_counts())
                
                if 'DKCC' in adata.obs:
                    print(f"\nCell types identified: {adata.obs['DKCC'].nunique()}")
                    print(adata.obs['DKCC'].value_counts().head(10))
                    if adata.obs['DKCC'].nunique() > 10:
                        print(f"... and {adata.obs['DKCC'].nunique() - 10} more")
                
                print("\n" + "="*60)
                print("[OK] Classification complete!")
                print("="*60)
        
        finally:
            # Cleanup temp files
            if self.verbose:
                print(f"\nCleaning up temporary files...")
            
            import shutil
            try:
                shutil.rmtree(temp_dir)
                if self.verbose:
                    print("  [OK] Cleanup complete")
            except Exception as e:
                if self.verbose:
                    print(f"  Warning: Could not remove temp directory: {e}")
        
        return adata
    
    def get_marker_genes(self, cell_type: str) -> pd.DataFrame:
        """
        Get marker genes for a specific cell type from the R package.
        
        Parameters
        ----------
        cell_type : str
            Cell type name
            
        Returns
        -------
        pd.DataFrame
            Marker genes for the specified cell type
        """
        if not self._r_packages_loaded:
            raise RuntimeError("R packages not loaded.")
        
        # This is a placeholder - implement based on what DevKidCC exposes
        warnings.warn("get_marker_genes not yet implemented")
        return pd.DataFrame()


def classify_kidney_cells(adata: ad.AnnData, 
                         copy: bool = True, 
                         verbose: bool = True,
                         install_deps: bool = True) -> ad.AnnData:
    """
    Classify kidney cells using DevKidCC (convenience wrapper).
    
    This function provides a simple one-line interface to classify
    kidney cells. It handles all R package setup automatically.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Single-cell RNA-seq data (AnnData format)
    copy : bool, default=True
        Whether to return a copy of the data
    verbose : bool, default=True
        Whether to print progress messages
    install_deps : bool, default=True
        Whether to automatically install missing R packages
        
    Returns
    -------
    anndata.AnnData
        Data with cell type classifications in .obs
        
    Examples
    --------
    >>> import scanpy as sc
    >>> from devkidcc import classify_kidney_cells
    >>> 
    >>> # Load your data
    >>> adata = sc.read_h5ad("kidney_organoid.h5ad")
    >>> 
    >>> # Classify (one line!)
    >>> adata = classify_kidney_cells(adata)
    >>> 
    >>> # View results
    >>> print(adata.obs['DKCC'].value_counts())
    """
    classifier = DevKidCCClassifier(verbose=verbose, install_deps=install_deps)
    return classifier.classify(adata, copy=copy)