"""Library sizes must be summed over the whole matrix, not the projection.

Seurat's LogNormalize divides each cell by its total counts, and that total is
summed over whatever genes are in the object. Because the wrapper projects onto
the DevKidCC reference genes before handing anything to R, letting R compute the
totals itself would divide by a *subset* sum -- shifting every normalised value
and every scPred score with it.

Measured on Howden 2019 (5,365 cells) against labels produced before the
projection existed: LineageID kappa 0.66 without this fix, 0.9989 with it. So
the column below is not bookkeeping; it is the difference between reproducing
the published labels and quietly not.
"""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from devkidcc.classifier_subprocess import _LIB_SIZE_COL, DevKidCCClassifier


@pytest.fixture
def classifier():
    obj = DevKidCCClassifier.__new__(DevKidCCClassifier)
    obj.verbose = False
    return obj


@pytest.fixture
def counts():
    """4 cells x 6 genes, only three of which are reference genes."""
    return np.array(
        [[0, 5, 1, 0, 2, 9],
         [7, 0, 0, 3, 1, 0],
         [2, 2, 2, 2, 2, 2],
         [0, 0, 0, 0, 0, 4]],
        dtype=np.float32,
    )


def _adata(counts, sparse):
    return ad.AnnData(
        X=csr_matrix(counts) if sparse else counts,
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(counts.shape[0])]),
        var=pd.DataFrame(index=["PAX2", "ZZZFAKE1", "SIX1", "ZZZFAKE2",
                                "GATA3", "ZZZFAKE3"]),
    )


@pytest.mark.parametrize("sparse", [True, False], ids=["sparse", "dense"])
def test_totals_cover_every_gene_not_just_reference_genes(classifier, counts, sparse):
    totals = classifier._library_sizes(_adata(counts, sparse))

    np.testing.assert_allclose(totals, counts.sum(axis=1))

    # The point of the test: the reference-gene subset sums to something else.
    reference_only = counts[:, [0, 2, 4]].sum(axis=1)
    assert not np.allclose(totals, reference_only)


def test_row_chunking_does_not_change_the_answer(classifier, monkeypatch):
    """A chunk size that does not divide the cell count evenly still sums right."""
    monkeypatch.setattr("devkidcc.classifier_subprocess._LIB_SIZE_ROW_CHUNK", 3)
    rng = np.random.default_rng(0)
    counts = rng.poisson(2.0, size=(11, 6)).astype(np.float32)

    totals = classifier._library_sizes(_adata(counts, sparse=False))

    np.testing.assert_allclose(totals, counts.sum(axis=1))


def test_r_script_normalises_from_the_column():
    """run_dkcc.R must read the column, not fall through to NormalizeData()."""
    from pathlib import Path

    import devkidcc

    source = (Path(devkidcc.__file__).parent / "run_dkcc.R").read_text()
    assert _LIB_SIZE_COL in source
    # The hand-rolled LogNormalize, and the stock call kept as the fallback for
    # callers driving the script with an unprojected matrix.
    assert "log1p" in source
    assert "NormalizeData(seurat_obj)" in source
