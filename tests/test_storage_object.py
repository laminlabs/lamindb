import pytest

import lamindb as ln
from lamindb.dev.storage.object import infer_suffix, write_to_file


def test_infer_suffix():
    import anndata as ad

    adata = ad.AnnData()
    with pytest.raises(ValueError):
        infer_suffix(adata, adata_format="my format")
    with pytest.raises(NotImplementedError):
        infer_suffix(ln.File)


def test_write_to_file():
    with pytest.raises(NotImplementedError):
        write_to_file(ln.File, "path")
