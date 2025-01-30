from pathlib import Path

import anndata as ad
import lamindb as ln
import pytest


@pytest.fixture(scope="module")
def local_filepath():
    return ln.core.datasets.anndata_file_pbmc68k_test().resolve()


@pytest.fixture(scope="module")
def adata(local_filepath):
    return ad.read_h5ad(local_filepath)


@pytest.fixture(scope="module")
def html_filepath():
    filepath = Path("./tmp.html")
    with open(filepath, "w") as f:
        f.write("<html><body><h1>Test</h1></body></html>")
    yield filepath
    filepath.unlink()


@pytest.fixture(scope="module")
def json_filepath():
    filepath = Path("./tmp.json")
    with open(filepath, "w") as f:
        f.write('{"a": 1}')
    yield filepath
    filepath.unlink()


def test_load_anndata(local_filepath, adata):
    artifact = ln.Artifact(local_filepath, description="test")
    assert local_filepath == artifact._local_filepath
    assert local_filepath == artifact.path
    assert local_filepath == artifact.cache()

    artifact = ln.Artifact.from_anndata(adata, description="test")
    assert artifact._memory_rep is adata
    assert artifact.load() is adata
    assert artifact._local_filepath.resolve() == artifact.cache() == artifact.path


def test_load_html(html_filepath):
    artifact = ln.Artifact(html_filepath, key=str(html_filepath))
    artifact.load()


def test_load_json(json_filepath):
    artifact = ln.Artifact(json_filepath, key=str(json_filepath))
    dictionary = artifact.load()
    assert dictionary["a"] == 1
