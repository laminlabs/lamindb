import json
import pickle
import shutil
from pathlib import Path

import anndata as ad
import lamindb as ln
import numpy as np
import pandas as pd
import pytest
from lamindb.core.storage._zarr import identify_zarr_type
from lamindb_setup.core.upath import (
    CloudPath,
)


@pytest.fixture(scope="session")
def get_small_adata():
    return ad.AnnData(
        X=np.array([[1, 2, 3], [4, 5, 6]]),
        obs={"feat1": ["A", "B"]},
        var=pd.DataFrame(index=["MYC", "TCF7", "GATA1"]),
        obsm={"X_pca": np.array([[1, 2], [3, 4]])},
    )


def test_zarr_upload_cache(get_small_adata):
    previous_storage = ln.setup.settings.storage.root_as_str
    ln.settings.storage = "s3://lamindb-test/core"

    zarr_path = Path("./test_adata.zarr")
    get_small_adata.write_zarr(zarr_path)

    artifact = ln.Artifact(zarr_path, key="test_adata.zarr")
    assert not artifact._storage_ongoing
    assert artifact.otype == "AnnData"
    assert artifact.n_files >= 1
    artifact.save()

    assert ln.Artifact.get(path=artifact.path) == artifact

    assert not artifact._storage_ongoing

    assert isinstance(artifact.path, CloudPath)
    assert artifact.path.exists()
    assert identify_zarr_type(artifact.path) == "anndata"

    shutil.rmtree(artifact.cache())

    cache_path = artifact._cache_path
    assert isinstance(artifact.load(), ad.AnnData)
    assert cache_path.is_dir()

    shutil.rmtree(cache_path)
    assert not cache_path.exists()
    artifact.cache()
    assert cache_path.is_dir()

    artifact.delete(permanent=True, storage=True)
    shutil.rmtree(zarr_path)

    # test zarr from memory
    artifact = ln.Artifact(get_small_adata, key="test_adata.anndata.zarr")
    assert not artifact._storage_ongoing
    assert artifact._local_filepath.is_dir()
    assert artifact.otype == "AnnData"
    assert artifact.suffix == ".anndata.zarr"
    assert artifact.n_files >= 1

    ln.save([artifact])  # use bulk save here for testing
    assert not artifact._storage_ongoing
    assert isinstance(artifact.path, CloudPath)
    assert artifact.path.exists()
    cache_path = artifact._cache_path
    assert cache_path.is_dir()

    shutil.rmtree(cache_path)
    assert not cache_path.exists()

    artifact._memory_rep = None

    assert isinstance(artifact.load(), ad.AnnData)
    assert cache_path.is_dir()

    artifact.delete(permanent=True, storage=True)

    ln.settings.storage = previous_storage


@pytest.mark.parametrize(
    "codec_version, allowed",
    [
        ("0.15.1", False),
        ("0.15.2", True),
        ("0.16.5", False),
        ("0.16.6", True),
        ("0.17.0", False),
        ("0.17.1", True),
    ],
)
def test_zarr_requires_patched_numcodecs(monkeypatch, tmp_path, codec_version, allowed):
    from lamindb.core.storage import _zarr

    monkeypatch.setattr(_zarr, "get_version", lambda package: codec_version)
    if allowed:
        assert _zarr.get_zarr_store(tmp_path) == str(tmp_path)
    else:
        with pytest.raises(ImportError, match="security-patched numcodecs"):
            _zarr.get_zarr_store(tmp_path)


@pytest.mark.parametrize("backed", [False, True])
def test_zarr_rejects_pickle(tmp_path, get_small_adata, backed):
    import numcodecs.errors
    import zarr
    from lamindb.core.loaders import load_to_memory
    from lamindb.core.storage._backed_access import backed_access
    from lamindb.core.storage._zarr import IS_ZARR_V3
    from lamindb_setup.core.upath import UPath

    path = tmp_path / "data.zarr"
    if IS_ZARR_V3:
        with ad.settings.override(zarr_write_format=2):
            get_small_adata.write_zarr(path)
    else:
        get_small_adata.write_zarr(path)

    def read():
        if backed:
            with backed_access(UPath(path)) as accessor:
                return accessor.X, accessor.uns
        adata = load_to_memory(path)
        return adata.X, adata.uns

    # Ordinary numeric and string data still read through the same entry points.
    np.testing.assert_array_equal(read()[0], get_small_adata.X)
    marker = tmp_path / "executed"

    class Payload:
        def __reduce__(self):
            return Path.touch, (marker,)

    # A v2 object array whose codec and chunk are both controlled by the store.
    evil = path / "uns" / "evil"
    evil.mkdir(parents=True)
    (evil / ".zarray").write_text(
        json.dumps(
            {
                "zarr_format": 2,
                "shape": [1],
                "chunks": [1],
                "dtype": "|O",
                "compressor": None,
                "fill_value": None,
                "order": "C",
                "filters": [{"id": "pickle", "protocol": 4}],
            }
        )
    )
    (evil / ".zattrs").write_text(
        json.dumps(
            {
                "encoding-type": "array",
                "encoding-version": "0.2.0",
            }
        )
    )
    (evil / "0").write_bytes(
        pickle.dumps(np.array([Payload()], dtype=object), protocol=4)
    )
    # Force readers to see the new array instead of older consolidated metadata.
    (path / ".zmetadata").unlink(missing_ok=True)
    disabled = getattr(
        zarr.errors, "PickleCodecDisabledError", numcodecs.errors.UnsafeDecodeError
    )
    # Zarr 3 can reject pickle object dtypes while parsing metadata, before decode.
    errors = (numcodecs.errors.UnsafeDecodeError, disabled, ValueError)
    with pytest.raises(errors, match="[Pp]ickle"):
        read()
    assert not marker.exists()
