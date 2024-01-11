from pathlib import Path

import lamindb as ln


def test_load_anndata():
    import anndata as ad

    local_filepath = ln.dev.datasets.anndata_file_pbmc68k_test().resolve()

    artifact = ln.Artifact(local_filepath, description="test")
    assert local_filepath == artifact._local_filepath
    assert local_filepath == artifact.path
    assert local_filepath == artifact.stage()

    adata = ad.read(local_filepath)
    artifact = ln.Artifact(adata, description="test")
    assert artifact._memory_rep is adata
    assert artifact.load() is adata
    assert artifact._local_filepath.resolve() == artifact.stage() == artifact.path


def test_load_html():
    filepath = Path("./tmp.html")
    with open(filepath, "w") as f:
        f.write("<html><body><h1>Test</h1></body></html>")
    artifact = ln.Artifact(filepath, key=str(filepath))
    artifact.load()
    filepath.unlink()


def test_load_json():
    filepath = Path("./tmp.json")
    with open(filepath, "w") as f:
        f.write('{"a": 1}')
    artifact = ln.Artifact(filepath, key=str(filepath))
    dictionary = artifact.load()
    assert dictionary["a"] == 1
    filepath.unlink()
