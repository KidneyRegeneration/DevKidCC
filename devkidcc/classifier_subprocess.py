"""
DevKidCC Python Wrapper (Subprocess-based)
Provides a Python interface to the R DevKidCC package using subprocess
This avoids the rpy2/reticulate segfault issue
"""

import warnings
import tempfile
import os
import subprocess
import shutil
from pathlib import Path
from typing import Optional, Union
import numpy as np
import pandas as pd
import anndata as ad


class DevKidCCClassifier:
    """
    Python wrapper for the DevKidCC R package (subprocess-based).

    This class provides a Python interface to classify kidney cells using
    the R DevKidCC package. It uses subprocess to call R, avoiding the
    rpy2/reticulate conflict that causes segfaults.

    Parameters
    ----------
    verbose : bool, default=True
        Whether to print progress messages
    rscript_path : str, optional
        Path to Rscript executable. If None, assumes 'Rscript' is in PATH

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

    def __init__(self, verbose: bool = True, rscript_path: Optional[str] = None):
        self.verbose = verbose
        self.rscript_path = rscript_path or 'Rscript'

        if self.verbose:
            print("Initializing DevKidCC Python wrapper (subprocess mode)...")

        # Find the R script
        script_dir = Path(__file__).parent
        self.r_script = script_dir / "run_dkcc.R"

        if not self.r_script.exists():
            raise FileNotFoundError(
                f"R script not found: {self.r_script}\n"
                "Please ensure run_dkcc.R is in the same directory as classifier_subprocess.py"
            )

        # Test that Rscript is available
        self._check_rscript()

        # Check R package dependencies
        self._check_r_packages()

        if self.verbose:
            print("[OK] Classifier initialized\n")

    def _check_rscript(self):
        """Check that Rscript is available."""
        try:
            result = subprocess.run(
                [self.rscript_path, '--version'],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode != 0:
                raise RuntimeError(f"Rscript not working: {result.stderr}")

            if self.verbose:
                # Extract version from output
                version_line = result.stdout.strip().split('\n')[0]
                print(f"  Found R: {version_line}")

        except FileNotFoundError:
            raise RuntimeError(
                f"Rscript not found at: {self.rscript_path}\n"
                "Please install R or provide the path to Rscript via rscript_path parameter"
            )
        except subprocess.TimeoutExpired:
            raise RuntimeError("Rscript --version timed out")

    def _check_r_packages(self):
        """Check that required R packages are installed."""
        if self.verbose:
            print("\n  Checking R package dependencies...")

        required_packages = ['Seurat', 'SeuratDisk', 'DevKidCC']

        for pkg in required_packages:
            check_code = f'if (!requireNamespace("{pkg}", quietly = TRUE)) quit(status = 1)'

            try:
                result = subprocess.run(
                    [self.rscript_path, '-e', check_code],
                    capture_output=True,
                    text=True,
                    timeout=30
                )

                if result.returncode != 0:
                    raise RuntimeError(
                        f"Required R package '{pkg}' not found.\n"
                        f"Please install it in R:\n"
                        f"  install.packages('{pkg}')  # For CRAN packages\n"
                        f"  remotes::install_github('...')  # For GitHub packages"
                    )

                if self.verbose:
                    print(f"    [OK] {pkg}")

            except subprocess.TimeoutExpired:
                raise RuntimeError(f"Timeout while checking for R package: {pkg}")

        if self.verbose:
            print("  [OK] All R packages found")

    def classify(self, adata: ad.AnnData, copy: bool = True) -> ad.AnnData:
        """
        Classify kidney cells using DevKidCC.

        This is the main classification method. It saves the AnnData object
        to a temporary h5ad file, calls the R script via subprocess to run
        DevKidCC classification, and reads the results back.

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
        if copy:
            adata = adata.copy()

        if self.verbose:
            print("=" * 60)
            print("DevKidCC Classification Pipeline (Subprocess Backend)")
            print("=" * 60)
            print(f"Input: {adata.n_obs} cells x {adata.n_vars} genes\n")

        # Create temp directory for files
        temp_dir = tempfile.mkdtemp()

        try:
            # Save expression matrix and metadata to CSV
            input_csv = os.path.join(temp_dir, "counts.csv")
            obs_csv = os.path.join(temp_dir, "obs.csv")
            output_csv = os.path.join(temp_dir, "results.csv")

            if self.verbose:
                print("Saving data to temporary files...")

            # Save count matrix (genes x cells, transposed for R)
            count_df = pd.DataFrame(
                adata.X.T.toarray() if hasattr(adata.X, 'toarray') else adata.X.T,
                index=adata.var_names,
                columns=adata.obs_names
            )
            count_df.to_csv(input_csv)

            # Save obs metadata
            adata.obs.to_csv(obs_csv)

            if self.verbose:
                print("  [OK] Data saved\n")

            # Call R script
            if self.verbose:
                print("Running DevKidCC classification in R...")
                print("(Output from R script will appear below)\n")
                print("-" * 60)

            result = subprocess.run(
                [self.rscript_path, str(self.r_script), input_csv, output_csv, obs_csv],
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout
            )

            if self.verbose:
                print(result.stdout)
                print("-" * 60)

            if result.returncode != 0:
                error_msg = f"R script failed with exit code {result.returncode}\n"
                if result.stderr:
                    error_msg += f"Error output:\n{result.stderr}"
                raise RuntimeError(error_msg)

            # Load results
            if not os.path.exists(output_csv):
                raise RuntimeError(
                    f"Output file not created: {output_csv}\n"
                    "R script may have failed silently"
                )

            if self.verbose:
                print("\nLoading results back to Python...")

            # Read results CSV
            result_metadata = pd.read_csv(output_csv, index_col=0)

            # Ensure index matches
            if not result_metadata.index.equals(adata.obs.index):
                # Try to align by cell names
                result_metadata = result_metadata.loc[adata.obs.index]

            # Copy all columns from results to adata.obs
            for col in result_metadata.columns:
                adata.obs[col] = result_metadata[col]

            if self.verbose:
                print("  [OK] Results loaded")

                print("\n" + "=" * 60)
                print("Classification Summary")
                print("=" * 60)

                if 'LineageID' in adata.obs:
                    print(f"\nLineages identified: {adata.obs['LineageID'].nunique()}")
                    print(adata.obs['LineageID'].value_counts())

                if 'DKCC' in adata.obs:
                    print(f"\nCell types identified: {adata.obs['DKCC'].nunique()}")
                    print("Top 10 cell types:")
                    print(adata.obs['DKCC'].value_counts().head(10))
                    if adata.obs['DKCC'].nunique() > 10:
                        print(f"... and {adata.obs['DKCC'].nunique() - 10} more")

                print("\n" + "=" * 60)
                print("[OK] Classification complete!")
                print("=" * 60)

        except subprocess.TimeoutExpired:
            raise RuntimeError(
                "R script timed out after 1 hour.\n"
                "This may happen with very large datasets."
            )

        finally:
            # Cleanup temp files
            if self.verbose:
                print("\nCleaning up temporary files...")

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
        warnings.warn("get_marker_genes not yet implemented in subprocess mode")
        return pd.DataFrame()


def classify_kidney_cells(adata: ad.AnnData,
                         copy: bool = True,
                         verbose: bool = True,
                         rscript_path: Optional[str] = None) -> ad.AnnData:
    """
    Classify kidney cells using DevKidCC (convenience wrapper).

    This function provides a simple one-line interface to classify
    kidney cells using the subprocess-based wrapper.

    Parameters
    ----------
    adata : anndata.AnnData
        Single-cell RNA-seq data (AnnData format)
    copy : bool, default=True
        Whether to return a copy of the data
    verbose : bool, default=True
        Whether to print progress messages
    rscript_path : str, optional
        Path to Rscript executable

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
    classifier = DevKidCCClassifier(verbose=verbose, rscript_path=rscript_path)
    return classifier.classify(adata, copy=copy)
