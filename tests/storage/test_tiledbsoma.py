import shutil
from pathlib import Path

import lamindb as ln
import numpy as np
import pytest
import tiledbsoma
import tiledbsoma.io
from lamindb.core.loaders import load_h5ad
from lamindb.core.storage._tiledbsoma import (
    _open_tiledbsoma,
    _soma_store_n_observations,
    _tiledb_config_s3,
)
from lamindb.integrations import save_tiledbsoma_experiment


@pytest.mark.parametrize("storage", [None, "s3://lamindb-test/storage"])
def test_write_read_tiledbsoma(storage):
    if storage is not None:
        previous_storage = ln.setup.settings.storage.root_as_str
        ln.settings.storage = storage

    test_file = ln.examples.datasets.anndata_file_pbmc68k_test()
    adata = load_h5ad(test_file)
    # write less
    adata = adata[:5, :2].copy()
    del adata.varp
    del adata.obsp
    del adata.layers
    del adata.uns  # seems to cause problems for append
    if storage is None:
        # test local with zarr
        test_file = test_file.with_suffix(".zarr")
        adata.write_zarr(test_file)
    else:
        adata.write_h5ad(test_file)

    create_transform = ln.Transform(key="test create tiledbsoma store").save()
    create_run = ln.Run(create_transform).save()

    # fails with a view
    with pytest.raises(ValueError, match="Can not write an `AnnData` view"):
        save_tiledbsoma_experiment([adata[:2]], run=create_run, measurement_name="RNA")

    artifact_soma = save_tiledbsoma_experiment(
        [test_file],
        description="test tiledbsoma",
        key="scrna/my-big-dataset.tiledbsoma",  # can also be None, but that's trivial
        run=create_run,
        measurement_name="RNA",
    )
    assert artifact_soma.path.stem == artifact_soma.uid[:16]
    assert artifact_soma.key == "scrna/my-big-dataset.tiledbsoma"
    assert artifact_soma.suffix == ".tiledbsoma"
    assert artifact_soma._key_is_virtual
    assert artifact_soma.otype == "tiledbsoma"
    assert artifact_soma.n_observations == adata.n_obs

    with artifact_soma.open() as store:  # mode="r" by default
        assert isinstance(store, tiledbsoma.Experiment)
        obs = store["obs"]
        n_obs = len(obs)
        assert n_obs == adata.n_obs
        assert "lamin_run_uid" in obs.schema.names
        run_ids = (
            obs.read(column_names=["lamin_run_uid"])
            .concat()
            .to_pandas()["lamin_run_uid"]
        )
        assert all(run_ids == create_run.uid)
        assert set(run_ids.cat.categories) == {create_run.uid}
        # test reading X
        ms_rna = store.ms["RNA"]
        n_vars = len(ms_rna.var)
        assert n_vars == adata.n_vars
        X = ms_rna["X"]["data"].read().coos((n_obs, n_vars)).concat().to_scipy()
        assert X.sum() == adata.X.sum()

    cache_path = artifact_soma.cache()
    hash_before_changes = artifact_soma.hash
    with artifact_soma.open(mode="w") as store:
        assert store.__class__.__name__ == "ExperimentTrack"
        tiledbsoma.io.add_matrix_to_collection(
            exp=store,
            measurement_name="RNA",
            collection_name="obsm",
            matrix_name="test_array",
            matrix_data=np.ones((n_obs, 2)),
        )
    assert artifact_soma.hash != hash_before_changes
    assert artifact_soma.uid.endswith("0001")
    if storage is not None:
        # cache should be ignored and deleted after the changes
        assert not cache_path.exists()
    else:
        assert artifact_soma.path == cache_path

    adata_to_append_1 = adata[:3].copy()
    adata_to_append_1.obs["obs_id"] = adata_to_append_1.obs.index.to_numpy() + "***"
    adata_to_append_1.var["var_id"] = adata_to_append_1.var.index
    adata_to_append_2 = adata[3:5].copy()
    adata_to_append_2.obs["obs_id"] = adata_to_append_2.obs.index.to_numpy() + "***"
    adata_to_append_2.var["var_id"] = adata_to_append_2.var.index
    adata_to_append_2.write_h5ad("adata_to_append_2.h5ad")

    append_transform = ln.Transform(key="test append tiledbsoma store").save()
    append_run = ln.Run(append_transform).save()

    # here run should be passed
    with pytest.raises(ValueError, match="Pass `run`"):
        save_tiledbsoma_experiment(
            [adata_to_append_1],
            revises=artifact_soma,
            run=None,
            measurement_name="RNA",
        )

    artifact_soma_append = save_tiledbsoma_experiment(
        [adata_to_append_1, "adata_to_append_2.h5ad"],
        revises=artifact_soma,
        run=append_run,
        measurement_name="RNA",
        append_obsm_varm=True,
    )
    assert artifact_soma_append.uid.endswith("0002")
    # below is inherited from "scrna/my-big-dataset.tiledbsoma"
    assert artifact_soma_append.key == "scrna/my-big-dataset.tiledbsoma"

    # wrong mode, should be either r or w for tiledbsoma
    with pytest.raises(ValueError):
        artifact_soma_append.open(mode="p")

    # test running without the context manager
    store = artifact_soma_append.open()
    n_obs_final = adata.n_obs + sum(
        adt.n_obs for adt in [adata_to_append_1, adata_to_append_2]
    )
    obs = store["obs"]
    assert len(obs) == n_obs_final == artifact_soma_append.n_observations
    run_ids = (
        obs.read(column_names=["lamin_run_uid"])
        .concat()
        .to_pandas()["lamin_run_uid"]
        .cat.categories
    )
    assert set(run_ids) == {create_run.uid, append_run.uid}
    store.close()

    # test correctness of deletion for _overwrite_versions=True
    soma_path = artifact_soma_append.path
    assert soma_path.exists()
    # select specific version and delete
    # check that the store is stil there
    assert soma_path.exists()
    assert ln.Artifact.filter(description="test tiledbsoma").count() == 3
    artifact_soma_append.versions.filter(uid__endswith="0001").one().delete(
        permanent=True
    )
    assert soma_path.exists()
    assert ln.Artifact.filter(description="test tiledbsoma").count() == 2
    # make sure it the store is actually deleted
    artifact_soma_append.delete(permanent=True)
    assert not soma_path.exists()
    assert not ln.Artifact.filter(description="test tiledbsoma").exists()

    Path("adata_to_append_2.h5ad").unlink()

    if storage is not None:
        ln.settings.storage = previous_storage


def test_from_tiledbsoma():
    test_file = ln.examples.datasets.anndata_file_pbmc68k_test()
    soma_path = "mystore.tiledbsoma"
    tiledbsoma.io.from_h5ad(soma_path, test_file, measurement_name="RNA")
    # wrong suffix
    with pytest.raises(ValueError):
        ln.Artifact.from_tiledbsoma("mystore")

    artifact = ln.Artifact.from_tiledbsoma(
        soma_path, description="test soma store"
    ).save()
    assert artifact.n_observations == 30

    with _open_tiledbsoma(artifact.path, mode="r") as store:
        # experiment
        assert _soma_store_n_observations(store) == 30
        # dataframe
        assert _soma_store_n_observations(store.obs) == 30
        # treat as unstructured collection, data + raw
        assert _soma_store_n_observations(store.ms) == 60
        # measurement
        assert _soma_store_n_observations(store.ms["RNA"]) == 30
        # array
        assert _soma_store_n_observations(store.ms["RNA"]["X"]["data"]) == 30

    artifact.delete(permanent=True)
    shutil.rmtree(soma_path)


def test_tiledb_config():
    storepath = ln.UPath("s3://bucket/key?endpoint_url=http://localhost:9000/s3")
    tiledb_config = _tiledb_config_s3(storepath)
    assert tiledb_config["vfs.s3.endpoint_override"] == "localhost:9000/s3"
    assert tiledb_config["vfs.s3.scheme"] == "http"
    assert tiledb_config["vfs.s3.use_virtual_addressing"] == "false"
    assert tiledb_config["vfs.s3.region"] == ""
