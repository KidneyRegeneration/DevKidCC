"""End-to-end regression against a stored classification.

Opt-in and slow (tens of minutes): it runs R. Point it at an h5ad that already
carries reference labels and run with `-m slow`:

    DEVKIDCC_REGRESSION_H5AD=/path/Howden_2019_Organoids_classified.h5ad \\
    DEVKIDCC_REGRESSION_KNN_COL=DKCC_v05_knn \\
    DEVKIDCC_REGRESSION_NOKNN_COL=DKCC_v05_noknn \\
    TMPDIR=/somewhere/with/space \\
    pytest tests/test_regression_howden.py -m slow -s

The two arms are held to different bars deliberately. Without smoothing the
result should reproduce the stored labels almost exactly. With smoothing it may
drift a little: DKCC() now computes the UMAP the KNN vote runs over *itself*,
after the zero-variance filter, whereas the labels were produced by a script
that computed it beforehand on the unfiltered matrix. A few percent is that
change; a large gap is not, and means the old in-script override was doing
something the package is not.
"""

import os

import numpy as np
import pytest

pytestmark = pytest.mark.slow

KNN_MIN_AGREEMENT = 0.90
NOKNN_MIN_AGREEMENT = 0.99


def _load():
    path = os.environ.get("DEVKIDCC_REGRESSION_H5AD")
    if not path or not os.path.exists(path):
        pytest.skip("set DEVKIDCC_REGRESSION_H5AD to a classified h5ad")
    import anndata as ad
    return ad.read_h5ad(path)


def _agreement(observed, expected):
    observed = np.asarray(observed, dtype=object)
    expected = np.asarray(expected, dtype=object)
    return float((observed == expected).mean())


@pytest.mark.parametrize("arm,env_col,floor", [
    ("noknn", "DEVKIDCC_REGRESSION_NOKNN_COL", NOKNN_MIN_AGREEMENT),
    ("knn", "DEVKIDCC_REGRESSION_KNN_COL", KNN_MIN_AGREEMENT),
])
def test_matches_stored_labels(arm, env_col, floor):
    adata = _load()

    col = os.environ.get(env_col)
    if not col or col not in adata.obs:
        pytest.skip(f"${env_col} must name a column present in the h5ad")

    from devkidcc import classify_kidney_cells

    n_vars_before = adata.n_vars
    result = classify_kidney_cells(
        adata, knn_smooth=(arm == "knn"), copy=True, verbose=True
    )

    assert "DKCC" in result.obs
    assert "LineageID" in result.obs
    # The reference-gene projection is for the subprocess only.
    assert result.n_vars == n_vars_before

    agreement = _agreement(result.obs["DKCC"], adata.obs[col])
    print(f"\n[{arm}] DKCC agreement vs {col}: {agreement:.4f}")
    assert agreement >= floor, (
        f"{arm} arm agreed on only {agreement:.1%} of cells (floor {floor:.0%})"
    )
