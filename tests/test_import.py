"""Package surface: what devkidcc exports and what those callables accept.

These guard the API that downstream pipelines bind to. `classify()` gaining or
losing a keyword has silently broken callers before.
"""

import inspect

import pytest

import devkidcc


def test_exports():
    for name in ("DevKidCCClassifier", "classify_kidney_cells",
                 "load_reference_genes", "reference_gene_path"):
        assert name in devkidcc.__all__
        assert hasattr(devkidcc, name)


def test_no_rpy2_dependency():
    """The subprocess backend must not import rpy2 — that is its whole point.

    rpy2 and reticulate segfault when both are loaded into one process, which
    is why classification goes out over a subprocess instead.
    """
    import ast

    import devkidcc.classifier_subprocess as mod

    tree = ast.parse(inspect.getsource(mod))
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])

    assert not imported & {"rpy2", "anndata2ri"}


@pytest.mark.parametrize("expected", [
    "adata", "threshold", "max_iter", "copy", "knn_smooth", "knn_iter",
])
def test_classify_signature(expected):
    params = inspect.signature(devkidcc.DevKidCCClassifier.classify).parameters
    assert expected in params


@pytest.mark.parametrize("expected", [
    "adata", "threshold", "max_iter", "copy", "verbose", "rscript_path",
    "knn_smooth", "knn_iter",
])
def test_classify_kidney_cells_signature(expected):
    params = inspect.signature(devkidcc.classify_kidney_cells).parameters
    assert expected in params


def test_classify_defaults_smooth_on():
    params = inspect.signature(devkidcc.DevKidCCClassifier.classify).parameters
    assert params["knn_smooth"].default is True
    assert params["knn_iter"].default is None
    assert params["copy"].default is True


def test_constructor_takes_no_install_deps():
    """Historic tests passed install_deps=True; no constructor ever accepted it."""
    params = inspect.signature(devkidcc.DevKidCCClassifier.__init__).parameters
    assert "install_deps" not in params
    assert set(params) == {"self", "verbose", "rscript_path"}
