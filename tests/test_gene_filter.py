"""The reference-gene projection: it must shrink what R sees, and nothing else.

The bug these exist to prevent: the filter used to rebind `adata`, so callers
using copy=False got a truncated object back — ~10k genes in place of their
whole matrix, silently, in the middle of a pipeline.
"""

import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

import devkidcc
from devkidcc.classifier_subprocess import (
    DevKidCCClassifier,
    load_reference_genes,
    reference_gene_path,
)


@pytest.fixture
def classifier():
    """A classifier that skips __init__ — these tests never invoke R."""
    obj = DevKidCCClassifier.__new__(DevKidCCClassifier)
    obj.verbose = False
    return obj


@pytest.fixture
def adata():
    """4 cells x 6 genes; PAX2/SIX1/GATA3 are real reference genes, ZZZ* are not."""
    X = csr_matrix(np.arange(24, dtype=np.float32).reshape(4, 6))
    return ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(4)]),
        var=pd.DataFrame(index=["PAX2", "ZZZFAKE1", "SIX1", "ZZZFAKE2",
                                "GATA3", "ZZZFAKE3"]),
    )


def test_packaged_gene_list_is_present_and_plausible():
    genes = load_reference_genes()
    assert len(genes) > 9000
    assert {"PAX2", "SIX1", "GATA3"} <= genes


def test_missing_env_override_raises(monkeypatch, tmp_path):
    monkeypatch.setenv("DEVKIDCC_REF_GENES", str(tmp_path / "nope.txt"))
    with pytest.raises(FileNotFoundError):
        reference_gene_path()


def test_env_override_wins(monkeypatch, tmp_path):
    custom = tmp_path / "genes.txt"
    custom.write_text("PAX2\nSIX1\n")
    monkeypatch.setenv("DEVKIDCC_REF_GENES", str(custom))
    assert reference_gene_path() == custom
    assert load_reference_genes() == {"PAX2", "SIX1"}


def test_csv_written_holds_only_reference_genes(classifier, adata, tmp_path,
                                                monkeypatch):
    custom = tmp_path / "genes.txt"
    custom.write_text("PAX2\nSIX1\nGATA3\n")
    monkeypatch.setenv("DEVKIDCC_REF_GENES", str(custom))

    ref = load_reference_genes()
    keep = np.flatnonzero(np.array([g in ref for g in adata.var_names]))
    out = tmp_path / "counts.csv"
    classifier._write_counts_csv(adata, keep, str(out))

    written = pd.read_csv(out, index_col=0)
    assert list(written.index) == ["PAX2", "SIX1", "GATA3"]
    assert list(written.columns) == list(adata.obs_names)


def test_csv_values_are_the_transposed_matrix(classifier, adata, tmp_path):
    keep = np.arange(adata.n_vars)
    out = tmp_path / "counts.csv"
    classifier._write_counts_csv(adata, keep, str(out))

    written = pd.read_csv(out, index_col=0)
    np.testing.assert_allclose(written.to_numpy(), adata.X.toarray().T)


def test_dense_matrix_writes_identically(classifier, adata, tmp_path):
    dense = ad.AnnData(X=adata.X.toarray(), obs=adata.obs, var=adata.var)
    keep = np.arange(dense.n_vars)

    sparse_out, dense_out = tmp_path / "s.csv", tmp_path / "d.csv"
    classifier._write_counts_csv(adata, keep, str(sparse_out))
    classifier._write_counts_csv(dense, keep, str(dense_out))

    np.testing.assert_allclose(
        pd.read_csv(sparse_out, index_col=0).to_numpy(),
        pd.read_csv(dense_out, index_col=0).to_numpy(),
    )


def test_chunk_boundary_is_respected(classifier, tmp_path, monkeypatch):
    """More genes than one write batch — off-by-one here loses rows silently."""
    monkeypatch.setattr("devkidcc.classifier_subprocess._CSV_GENE_CHUNK", 7)
    n_genes = 23
    a = ad.AnnData(
        X=csr_matrix(np.arange(2 * n_genes, dtype=np.float32).reshape(2, n_genes)),
        obs=pd.DataFrame(index=["c0", "c1"]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]),
    )
    out = tmp_path / "counts.csv"
    classifier._write_counts_csv(a, np.arange(n_genes), str(out))

    written = pd.read_csv(out, index_col=0)
    assert list(written.index) == list(a.var_names)
    np.testing.assert_allclose(written.to_numpy(), a.X.toarray().T)


def test_write_does_not_mutate_or_truncate_the_caller_object(classifier, adata,
                                                             tmp_path):
    n_vars_before, n_obs_before = adata.n_vars, adata.n_obs
    var_before = list(adata.var_names)

    classifier._write_counts_csv(adata, np.array([0, 2]), str(tmp_path / "c.csv"))

    assert adata.n_vars == n_vars_before
    assert adata.n_obs == n_obs_before
    assert list(adata.var_names) == var_before
