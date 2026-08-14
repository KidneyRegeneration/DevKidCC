"""
DevKidCC Python Wrapper (Subprocess-based)
Provides a Python interface to the R DevKidCC package using subprocess
This avoids the rpy2/reticulate segfault issue
"""

import gc
import warnings
import tempfile
import os
import subprocess
import shutil
import sys
from pathlib import Path
from typing import Optional, Union
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import issparse, csc_matrix

# Genes per batch when streaming the counts matrix out to CSV.
_CSV_GENE_CHUNK = 500

# Iterations passed to DKCC()'s knn.iter when KNN smoothing is on. Matches the
# R-side default; 0 disables smoothing entirely.
_DEFAULT_KNN_ITER = 20

# Cells per batch when summing library sizes from a dense matrix.
_LIB_SIZE_ROW_CHUNK = 5000

# Reserved obs column carrying each cell's total counts over the *full* input
# matrix. run_dkcc.R normalises with these instead of the projected column sums,
# so the reference-gene projection cannot move the scPred scores. Named to be
# unmistakably ours; it is stripped from the returned object.
_LIB_SIZE_COL = "dkcc_full_library_size"


def reference_gene_path() -> Path:
    """
    Locate the DevKidCC reference gene list.

    Resolution order: ``$DEVKIDCC_REF_GENES``, then the copy shipped inside the
    package. Raises if neither is present — falling back to the full matrix
    would silently reintroduce the out-of-memory kill the filter exists to
    prevent.
    """
    override = os.environ.get('DEVKIDCC_REF_GENES')
    if override:
        path = Path(override)
        if not path.exists():
            raise FileNotFoundError(
                f"$DEVKIDCC_REF_GENES points at {path}, which does not exist."
            )
        return path

    packaged = Path(__file__).parent / "data" / "reference_genes.txt"
    if packaged.exists():
        return packaged

    raise FileNotFoundError(
        "DevKidCC reference gene list not found. Expected it at "
        f"{packaged}, or set $DEVKIDCC_REF_GENES to a gene-per-line file. "
        "Regenerate it with scripts/export_reference_genes.R."
    )


def load_reference_genes() -> set:
    """Read the reference gene list as a set of symbols."""
    with open(reference_gene_path()) as fh:
        return {line.strip() for line in fh if line.strip()}


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

        required_packages = ['Seurat', 'DevKidCC']

        for pkg in required_packages:
            check_code = f'if (!requireNamespace("{pkg}", quietly = TRUE)) quit(status = 1)'

            try:
                result = subprocess.run(
                    [self.rscript_path, '-e', check_code],
                    capture_output=True,
                    text=True,
                    timeout=120  # Increased to handle renv dependency discovery (~40s on Windows)
                    # Note: First-time package checks may take longer due to R package loading
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

    def classify(self,
                 adata: ad.AnnData,
                 threshold: float = 0.7,
                 max_iter: int = 1,
                 copy: bool = True,
                 knn_smooth: bool = True,
                 knn_iter: Optional[int] = None) -> ad.AnnData:
        """
        Classify kidney cells using DevKidCC.

        This is the main classification method. It saves the AnnData object
        to a temporary h5ad file, calls the R script via subprocess to run
        DevKidCC classification, and reads the results back.

        Parameters
        ----------
        adata : anndata.AnnData
            Input single-cell RNA-seq data
        threshold : float, default=0.7
            Confidence threshold for cell type assignment (0-1).
            Higher values = stricter classification (fewer cells assigned).
            Lower values = more lenient (more cells assigned).
        max_iter : int, default=1
            Maximum number of iterations for refinement.
            Usually 1 is sufficient.
        copy : bool, default=True
            Whether to return a copy (recommended)
        knn_smooth : bool, default=True
            Whether to rescue unassigned cells by KNN vote over the UMAP
            embedding. Set False to benchmark the raw scPred assignments.
        knn_iter : int, optional
            Explicit iteration count for the KNN rescue, overriding
            `knn_smooth`. Maps directly onto DKCC()'s `knn.iter`; 0 disables.

        Returns
        -------
        anndata.AnnData
            Annotated data with classifications in .obs:
            - 'LineageID': Broad lineage category
            - 'DKCC': Detailed cell type annotation

            The returned object keeps every gene it came in with — the
            reference-gene filter below applies only to what is handed to R.

        Examples
        --------
        >>> import scanpy as sc
        >>> from devkidcc import DevKidCCClassifier
        >>>
        >>> adata = sc.read_h5ad("kidney_data.h5ad")
        >>> classifier = DevKidCCClassifier()
        >>>
        >>> # Default parameters
        >>> adata = classifier.classify(adata)
        >>>
        >>> # Stricter classification
        >>> adata = classifier.classify(adata, threshold=0.9)
        >>>
        >>> # More lenient classification
        >>> adata = classifier.classify(adata, threshold=0.5)
        """
        if copy:
            adata = adata.copy()

        if knn_iter is None:
            knn_iter = _DEFAULT_KNN_ITER if knn_smooth else 0

        if self.verbose:
            print("=" * 60)
            print("DevKidCC Classification Pipeline (Subprocess Backend)")
            print("=" * 60)
            print(f"Input: {adata.n_obs} cells x {adata.n_vars} genes")
            print(f"Parameters: threshold={threshold}, max_iter={max_iter}, "
                  f"knn_iter={knn_iter}\n")

        # Create temp directory for files
        temp_dir = tempfile.mkdtemp()

        try:
            # Save expression matrix and metadata to CSV
            input_csv = os.path.join(temp_dir, "counts.csv")
            obs_csv = os.path.join(temp_dir, "obs.csv")
            output_csv = os.path.join(temp_dir, "results.csv")

            if self.verbose:
                print("Saving data to temporary files...")

            # Restrict what goes to R to the DevKidCC reference genes. The full
            # matrix (e.g. 56k x 40k) produces a ~12 GB CSV that R must load
            # entirely into RAM (~25 GB Seurat object) -> OOM kill. The ~9,977
            # reference genes shrink that to ~3 GB of CSV and a ~6 GB R peak.
            #
            # This is a projection for the subprocess only. `adata` is never
            # rebound: callers passing copy=False must get their own object back
            # with every gene still on it.
            ref_genes = load_reference_genes()
            keep_idx = np.flatnonzero(
                np.fromiter((g in ref_genes for g in adata.var_names),
                            dtype=bool, count=adata.n_vars)
            )
            if keep_idx.size == 0:
                raise ValueError(
                    "None of the input genes are DevKidCC reference genes. "
                    "Check that var_names are HGNC symbols and carry no genome "
                    "prefix (e.g. 'GRCh38_GAPDH')."
                )
            if self.verbose:
                print(f"  Gene filter: {keep_idx.size} / {adata.n_vars} "
                      "reference genes retained")

            self._write_counts_csv(adata, keep_idx, input_csv)

            # Seurat's LogNormalize divides each cell by its library size, and
            # that library size is summed over whatever genes are in the object.
            # Projecting onto the reference genes above therefore shifts every
            # normalised value, and with it every scPred score -- silently, and
            # by enough to matter (Howden 2019: LineageID kappa 0.66 against the
            # unprojected labels). Send the totals computed over the *full*
            # matrix so R can normalise as if no projection had happened.
            obs_out = adata.obs.copy()
            obs_out[_LIB_SIZE_COL] = self._library_sizes(adata)
            obs_out.to_csv(obs_csv)

            if self.verbose:
                print("  [OK] Data saved\n")

            # Call R script
            if self.verbose:
                print("Running DevKidCC classification in R...")
                print("(Output from R script will appear below)\n")
                print("-" * 60)

            # R's own progress goes straight to our stdout rather than being
            # buffered until the end: these runs take tens of minutes and a
            # silent terminal is indistinguishable from a hang.
            result = subprocess.run(
                [self.rscript_path, str(self.r_script), input_csv, output_csv, obs_csv,
                 str(threshold), str(max_iter), str(knn_iter)],
                stdout=sys.stdout,
                stderr=subprocess.PIPE,
                text=True,
                timeout=3600  # 1 hour timeout
            )

            if self.verbose:
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

            # Copy all columns from results to adata.obs. The library-size column
            # is ours, not a result: it went over to R with the metadata and comes
            # back on it, so drop it rather than leaving it on the caller's object.
            for col in result_metadata.columns:
                if col == _LIB_SIZE_COL:
                    continue
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

    def _write_counts_csv(self, adata: ad.AnnData, keep_idx: np.ndarray,
                          input_csv: str) -> None:
        """
        Stream the counts matrix out as genes x cells CSV, restricted to
        `keep_idx` columns of `adata`.

        Written in gene batches rather than through pandas: a single
        `pd.DataFrame(adata.X.T.toarray()).to_csv()` materialises a dense array
        and a DataFrame of it at the same time, which is what makes large inputs
        fall over before R ever starts.
        """
        X_csc = csc_matrix(adata.X) if issparse(adata.X) else None
        var_names = adata.var_names

        with open(input_csv, 'w') as fh:
            fh.write(',' + ','.join(adata.obs_names) + '\n')
            for start in range(0, keep_idx.size, _CSV_GENE_CHUNK):
                cols = keep_idx[start:start + _CSV_GENE_CHUNK]
                if X_csc is not None:
                    block = X_csc[:, cols].T.toarray()   # (chunk x cells)
                else:
                    block = np.asarray(adata.X[:, cols]).T
                for i, col in enumerate(cols):
                    fh.write(var_names[col] + ',' +
                             ','.join(map(str, block[i].tolist())) + '\n')
                del block
                gc.collect()

        del X_csc
        gc.collect()

    @staticmethod
    def _library_sizes(adata: ad.AnnData) -> np.ndarray:
        """
        Per-cell total counts over every gene in `adata`, summed in row blocks so
        a dense input is never copied whole.
        """
        if issparse(adata.X):
            return np.asarray(adata.X.sum(axis=1)).ravel()

        totals = np.empty(adata.n_obs, dtype=np.float64)
        for start in range(0, adata.n_obs, _LIB_SIZE_ROW_CHUNK):
            stop = min(start + _LIB_SIZE_ROW_CHUNK, adata.n_obs)
            totals[start:stop] = np.asarray(adata.X[start:stop]).sum(axis=1)
        return totals

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
                         threshold: float = 0.7,
                         max_iter: int = 1,
                         copy: bool = True,
                         verbose: bool = True,
                         rscript_path: Optional[str] = None,
                         knn_smooth: bool = True,
                         knn_iter: Optional[int] = None) -> ad.AnnData:
    """
    Classify kidney cells using DevKidCC (convenience wrapper).

    This function provides a simple one-line interface to classify
    kidney cells using the subprocess-based wrapper.

    Parameters
    ----------
    adata : anndata.AnnData
        Single-cell RNA-seq data (AnnData format)
    threshold : float, default=0.7
        Confidence threshold for cell type assignment (0-1)
    max_iter : int, default=1
        Maximum number of iterations for refinement
    copy : bool, default=True
        Whether to return a copy of the data
    verbose : bool, default=True
        Whether to print progress messages
    rscript_path : str, optional
        Path to Rscript executable
    knn_smooth : bool, default=True
        Whether to rescue unassigned cells by KNN vote over the UMAP embedding
    knn_iter : int, optional
        Explicit KNN iteration count, overriding `knn_smooth`; 0 disables

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
    >>> # Classify with defaults
    >>> adata = classify_kidney_cells(adata)
    >>>
    >>> # Stricter classification
    >>> adata = classify_kidney_cells(adata, threshold=0.9)
    >>>
    >>> # View results
    >>> print(adata.obs['DKCC'].value_counts())
    """
    classifier = DevKidCCClassifier(verbose=verbose, rscript_path=rscript_path)
    return classifier.classify(adata, threshold=threshold, max_iter=max_iter,
                               copy=copy, knn_smooth=knn_smooth,
                               knn_iter=knn_iter)
