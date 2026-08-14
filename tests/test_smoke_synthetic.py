"""Slow smoke test: does the whole Python -> R -> Python round trip work?

Synthetic counts, so the labels mean nothing biologically. What it proves is
that R is reachable, DevKidCC is installed with a knn.iter-aware DKCC(), the CSV
handoff parses on both sides, and the results merge back onto .obs.

    pytest tests/test_smoke_synthetic.py -m slow -s
"""

import numpy as np
import pytest

pytestmark = pytest.mark.slow

N_CELLS = 300
N_GENES = 2000


@pytest.fixture
def synthetic():
    import anndata as ad
    import pandas as pd

    from devkidcc import load_reference_genes

    # Draw gene names from the real reference list: the wrapper projects onto
    # it before calling R, so "Gene_0"-style names would leave nothing to send.
    rng = np.random.default_rng(0)
    genes = sorted(load_reference_genes())
    assert len(genes) >= N_GENES
    chosen = list(rng.choice(genes, size=N_GENES, replace=False))

    X = rng.negative_binomial(5, 0.3, size=(N_CELLS, N_GENES)).astype(np.float32)
    return ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(N_CELLS)]),
        var=pd.DataFrame(index=chosen),
    )


@pytest.mark.parametrize("knn_smooth", [False, True])
def test_round_trip(synthetic, knn_smooth):
    from devkidcc import classify_kidney_cells

    result = classify_kidney_cells(synthetic, knn_smooth=knn_smooth,
                                   copy=True, verbose=True)

    assert "LineageID" in result.obs
    assert "DKCC" in result.obs
    assert result.n_obs == N_CELLS
    assert result.n_vars == N_GENES          # nothing lost to the gene filter
    assert result.obs["DKCC"].notna().any()


def test_copy_false_leaves_genes_alone(synthetic):
    """The historic bug: copy=False callers got back a gene-filtered object."""
    from devkidcc import classify_kidney_cells

    result = classify_kidney_cells(synthetic, knn_smooth=False, copy=False,
                                   verbose=False)

    assert result is synthetic
    assert synthetic.n_vars == N_GENES
    assert "DKCC" in synthetic.obs
