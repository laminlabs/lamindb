from inspect import signature

import anndata as ad
import bionty as bt
import lamindb as ln
import numpy as np
import pandas as pd
import pytest
from django.db.models.deletion import ProtectedError
from lamindb import _collection
from scipy.sparse import csr_matrix


@pytest.fixture(scope="module")
def df():
    return pd.DataFrame({"feat1": [1, 2], "feat2": [3, 4]})


@pytest.fixture(scope="module")
def adata():
    return ad.AnnData(
        X=np.array([[1, 2, 3], [4, 5, 6]]),
        obs={"feat1": ["A", "B"]},
        var=pd.DataFrame(index=["MYC", "TCF7", "GATA1"]),
        obsm={"X_pca": np.array([[1, 2], [3, 4]])},
    )


@pytest.fixture(scope="module")
def adata2():
    return ad.AnnData(
        X=np.array([[1, 2, 5], [4, 5, 8]]),
        obs={"feat1": ["A", "B"]},
        var=pd.DataFrame(index=["MYC", "TCF7", "GATA1"]),
        obsm={"X_pca": np.array([[1, 2], [3, 4]])},
    )


def test_signatures():
    # this seems currently the easiest and most transparent
    # way to test violations of the signature equality
    # the MockORM class is needed to get inspect.signature
    # to work
    class Mock:
        pass

    # class methods
    class_methods = []
    for name in class_methods:
        setattr(Mock, name, getattr(_collection, name))
        assert signature(getattr(Mock, name)) == _collection.SIGS.pop(name)
    # methods
    for name, sig in _collection.SIGS.items():
        assert signature(getattr(_collection, name)) == sig


def test_from_single_artifact(adata):
    bt.settings.organism = "human"
    features = ln.Feature.from_df(adata.obs)
    validated = ln.Feature.validate(
        [feature.name for feature in features], field="name"
    )
    ln.save([feature for (feature, valid) in zip(features, validated) if valid])
    artifact = ln.Artifact.from_anndata(adata, description="My adata")
    if not artifact._state.adding:
        artifact.delete(permanent=True)  # make sure we get a fresh one
        artifact = ln.Artifact.from_anndata(adata, description="My adata")
    with pytest.raises(ValueError) as error:
        ln.Collection(artifact, name="Test")
    assert str(error.exconly()).startswith(
        "ValueError: Not all artifacts are yet saved, please save them"
    )
    artifact.save()
    with pytest.raises(ValueError) as error:
        ln.Collection(artifact, artifact)
    assert str(error.exconly()).startswith(
        "ValueError: Only one non-keyword arg allowed: artifacts"
    )
    transform = ln.Transform(name="My test transform")
    transform.save()
    run = ln.Run(transform)
    run.save()
    collection = ln.Collection(artifact, name="My new collection", run=run)
    collection.save()
    assert collection.run.input_artifacts.get() == artifact
    collection.delete(permanent=True)
    artifact.delete(permanent=True)
    assert ln.Artifact.filter(id=artifact.id).one_or_none() is None


def test_edge_cases(df):
    with pytest.raises(ValueError) as error:
        ln.Collection(df, invalid_param=1)
    assert str(error.exconly()).startswith(
        "ValueError: Only artifacts, name, run, description, reference, reference_type, visibility can be passed, you passed: "
    )
    with pytest.raises(ValueError) as error:
        ln.Collection(1, name="Invalid")
    assert str(error.exconly()).startswith(
        "ValueError: Artifact or List[Artifact] is allowed."
    )
    artifact = ln.Artifact.from_df(df, description="Test artifact")
    assert artifact._state.adding
    with pytest.raises(ValueError) as error:
        ln.Collection([artifact])
    assert str(error.exconly()).startswith(
        "ValueError: Not all artifacts are yet saved, please save them"
    )
    artifact.save()
    with pytest.raises(ValueError) as error:
        ln.Collection([artifact, artifact])
    assert str(error.exconly()).startswith(
        "ValueError: Please pass artifacts with distinct hashes: these ones are"
        " non-unique"
    )
    artifact.delete(permanent=True, storage=True)


def test_from_inconsistent_artifacts(df, adata):
    artifact1 = ln.Artifact.from_df(df, description="My test")
    artifact1.save()
    artifact2 = ln.Artifact.from_anndata(adata, description="My test2")
    artifact2.save()
    collection = ln.Collection([artifact1, artifact2], name="Inconsistent")
    collection.save()
    # test idempotency of .save()
    collection.save()
    # create a run context
    ln.track(transform=ln.Transform(name="My test transform"))
    # can iterate over them
    collection.cache()
    assert set(ln.core.run_context.run.input_collections.all()) == {collection}
    # loading will throw an error here
    with pytest.raises(RuntimeError) as error:
        collection.load()
    assert str(error.exconly()).startswith(
        "RuntimeError: Can only load collections where all artifacts have the same suffix"
    )
    collection.delete(permanent=True)
    artifact1.delete(permanent=True)
    artifact2.delete(permanent=True)
    ln.core.run_context.run = None
    ln.core.run_context.transform = None


def test_from_consistent_artifacts(adata, adata2):
    artifact1 = ln.Artifact.from_anndata(adata, description="My test")
    artifact1.save()
    artifact2 = ln.Artifact.from_anndata(adata2, description="My test2")
    artifact2.save()
    transform = ln.Transform(name="My test transform")
    transform.save()
    run = ln.Run(transform)
    run.save()
    collection = ln.Collection([artifact1, artifact2], name="My test", run=run)
    assert collection._state.adding
    collection.save()
    assert set(collection.run.input_artifacts.all()) == {artifact1, artifact2}
    adata_joined = collection.load()
    assert "artifact_uid" in adata_joined.obs.columns
    assert artifact1.uid in adata_joined.obs.artifact_uid.cat.categories

    # re-run with hash-based lookup
    collection2 = ln.Collection([artifact1, artifact2], name="My test 1", run=run)
    assert not collection2._state.adding
    assert collection2.id == collection.id
    assert collection2.name == "My test 1"

    collection.delete(permanent=True)
    artifact1.delete(permanent=True)
    artifact2.delete(permanent=True)


def test_collection_mapped(adata, adata2):
    adata.strings_to_categoricals()
    adata.obs["feat2"] = adata.obs["feat1"]
    adata.layers["layer1"] = adata.X.copy()
    adata.layers["layer1"][0, 0] = 0
    artifact1 = ln.Artifact.from_anndata(adata, description="Part one")
    artifact1.save()
    adata2.X = csr_matrix(adata2.X)
    adata2.layers["layer1"] = adata2.X.copy()
    adata2.obs["feat2"] = adata2.obs["feat1"]
    artifact2 = ln.Artifact.from_anndata(adata2, description="Part two", format="zarr")
    artifact2.save()
    adata3 = adata2.copy()
    adata3.var_names = ["A", "B", "C"]
    artifact3 = ln.Artifact.from_anndata(adata3, description="Other vars")
    artifact3.save()
    collection = ln.Collection([artifact1, artifact2], name="Gather")
    collection.save()
    collection_outer = ln.Collection(
        [artifact1, artifact2, artifact3], name="Gather outer"
    )
    collection_outer.save()

    # test encoders
    with pytest.raises(ValueError):
        ls_ds = collection.mapped(encode_labels=["feat1"])
    with pytest.raises(ValueError):
        ls_ds = collection.mapped(obs_keys="feat1", encode_labels=["feat3"])
    with pytest.raises(ValueError):
        ls_ds = collection.mapped(obs_keys="feat1", unknown_label={"feat3": "Unknown"})
    with collection.mapped(obs_keys=["feat1", "feat2"], unknown_label="A") as ls_ds:
        assert ls_ds.encoders["feat1"]["A"] == -1
        assert ls_ds.encoders["feat1"]["B"] == 0
        assert ls_ds.encoders["feat2"]["A"] == -1
        assert ls_ds.encoders["feat2"]["B"] == 0
        assert ls_ds[0]["feat1"] == -1
        assert ls_ds[1]["feat1"] == 0
        assert ls_ds[0]["feat2"] == -1
        assert ls_ds[1]["feat2"] == 0
    with collection.mapped(
        obs_keys=["feat1", "feat2"], unknown_label={"feat1": "A"}
    ) as ls_ds:
        assert ls_ds.encoders["feat1"]["A"] == -1
        assert ls_ds.encoders["feat1"]["B"] == 0
        # can't predict order of elements in set
        A_enc = ls_ds.encoders["feat2"]["A"]
        B_enc = ls_ds.encoders["feat2"]["B"]
        assert A_enc in (0, 1)
        assert B_enc in (0, 1)
        assert A_enc != B_enc
        assert ls_ds[0]["feat1"] == -1
        assert ls_ds[1]["feat1"] == 0
        assert ls_ds[0]["feat2"] == A_enc
        assert ls_ds[1]["feat2"] == B_enc
    with collection.mapped(
        obs_keys=["feat1", "feat2"], unknown_label="A", encode_labels=["feat1"]
    ) as ls_ds:
        assert ls_ds.encoders["feat1"]["A"] == -1
        assert ls_ds.encoders["feat1"]["B"] == 0
        assert "feat2" not in ls_ds.encoders
        assert ls_ds[0]["feat1"] == -1
        assert ls_ds[1]["feat1"] == 0
        assert ls_ds[0]["feat2"] == "A"
        assert ls_ds[1]["feat2"] == "B"

    ls_ds = collection.mapped(obs_keys="feat1")
    assert not ls_ds.closed

    assert len(ls_ds) == 4
    assert len(ls_ds[0]) == 3 and len(ls_ds[2]) == 3
    assert len(ls_ds[0]["X"]) == 3
    assert np.array_equal(ls_ds[2]["X"], np.array([1, 2, 5]))
    weights = ls_ds.get_label_weights("feat1")
    assert all(weights[1:] == weights[0])
    weights = ls_ds.get_label_weights(["feat1", "feat2"])
    assert all(weights[1:] == weights[0])
    ls_ds.close()
    assert ls_ds.closed
    del ls_ds

    with collection.mapped(obs_keys="feat1", join="inner", dtype="float32") as ls_ds:
        assert not ls_ds.closed
        assert len(ls_ds) == 4
        assert len(ls_ds[0]) == 3 and len(ls_ds[2]) == 3
        assert str(ls_ds[0]["X"].dtype) == "float32"
        assert str(ls_ds[2]["X"].dtype) == "float32"
    assert ls_ds.closed

    ls_ds = collection.mapped(obs_keys="feat1", parallel=True)
    assert len(ls_ds[0]) == 3 and len(ls_ds[2]) == 3
    assert ls_ds[0]["_store_idx"] == 0
    assert ls_ds[2]["_store_idx"] == 1

    ls_ds = collection.mapped(
        layers_keys=["layer1"], obsm_keys=["X_pca"], obs_keys="feat1"
    )
    assert np.array_equal(ls_ds[0]["layer1"], np.array([0, 2, 3]))
    assert np.array_equal(ls_ds[2]["layer1"], np.array([1, 2, 5]))
    assert np.array_equal(ls_ds[2]["obsm_X_pca"], np.array([1, 2]))
    assert np.array_equal(ls_ds[3]["obsm_X_pca"], np.array([3, 4]))
    assert ls_ds.shape == (4, 3)
    assert ls_ds.original_shapes[0] == (2, 3) and ls_ds.original_shapes[1] == (2, 3)
    ls_ds.close()

    with collection.mapped(obs_keys="feat1", stream=True) as ls_ds:
        assert len(ls_ds[0]) == 3 and len(ls_ds[2]) == 3

    with pytest.raises(ValueError):
        with collection_outer.mapped(obs_keys="feat1", join="inner"):
            pass

    with collection_outer.mapped(
        layers_keys="X", obsm_keys="X_pca", obs_keys="feat1", join="outer"
    ) as ls_ds:
        assert ls_ds.shape == (6, 6)
        assert ls_ds.join_vars == "outer"
        assert len(ls_ds.var_joint) == 6
        assert len(ls_ds[0]) == 4
        assert len(ls_ds[0]["X"]) == 6
        assert np.array_equal(ls_ds[0]["X"], np.array([0, 0, 0, 3, 1, 2]))
        assert np.array_equal(ls_ds[1]["X"], np.array([0, 0, 0, 6, 4, 5]))
        assert np.array_equal(ls_ds[2]["X"], np.array([0, 0, 0, 5, 1, 2]))
        assert np.array_equal(ls_ds[3]["X"], np.array([0, 0, 0, 8, 4, 5]))
        assert np.array_equal(ls_ds[4]["X"], np.array([1, 2, 5, 0, 0, 0]))
        assert np.array_equal(ls_ds[5]["X"], np.array([4, 5, 8, 0, 0, 0]))
        assert np.issubdtype(ls_ds[2]["X"].dtype, np.integer)
        assert np.issubdtype(ls_ds[4]["X"].dtype, np.integer)
        assert np.array_equal(ls_ds[3]["obsm_X_pca"], np.array([3, 4]))

    with collection_outer.mapped(layers_keys="layer1", join="outer") as ls_ds:
        assert np.array_equal(ls_ds[0]["layer1"], np.array([0, 0, 0, 3, 0, 2]))
        assert np.array_equal(ls_ds[4]["layer1"], np.array([1, 2, 5, 0, 0, 0]))

    collection.delete(permanent=True)
    collection_outer.delete(permanent=True)
    artifact1.delete(permanent=True)
    artifact2.delete(permanent=True)
    artifact3.delete(permanent=True)


def test_is_new_version_of_versioned_collection(df, adata):
    # create a versioned collection
    artifact = ln.Artifact.from_df(df, description="test")
    artifact.save()
    collection = ln.Collection(artifact, name="test", version="1")
    assert collection.version == "1"
    collection.save()

    artifact = ln.Artifact.from_anndata(adata, description="test")
    artifact.save()

    with pytest.raises(ValueError) as error:
        collection_v2 = ln.Collection(
            artifact, is_new_version_of=collection, version="1"
        )
    assert error.exconly() == "ValueError: Please increment the previous version: '1'"

    # create new collection from old collection
    collection_v2 = ln.Collection(artifact, is_new_version_of=collection)
    assert collection.version == "1"
    assert collection_v2.stem_uid == collection.stem_uid
    assert collection_v2.version == "2"
    assert collection_v2.name == "test"

    collection_v2.save()

    # create new collection from newly versioned collection
    df.iloc[0, 0] = 0
    artifact = ln.Artifact.from_df(df, description="test")
    artifact.save()
    collection_v3 = ln.Collection(
        artifact, name="test1", is_new_version_of=collection_v2
    )
    assert collection_v3.stem_uid == collection.stem_uid
    assert collection_v3.version == "3"
    assert collection_v3.name == "test1"

    artifacts_v2 = collection_v2.artifacts.all()
    collection_v2.delete(permanent=True)
    artifacts_v2.delete(permanent=True)
    artifacts = collection.artifacts.all()
    collection.delete(permanent=True)
    artifacts.delete(permanent=True)


def test_is_new_version_of_unversioned_collection(df, adata):
    # unversioned collection
    artifact = ln.Artifact.from_df(df, description="test")
    artifact.save()
    collection = ln.Collection(artifact, name="test2")
    assert collection.version is None

    # what happens if we don't save the old collection?
    # add a test for it!
    collection.save()

    with pytest.raises(TypeError):
        ln.Collection(adata, is_new_version_of="wrong-type")

    artifact2 = ln.Artifact.from_anndata(adata, description="test")
    artifact2.save()

    # create new collection from old collection
    new_collection = ln.Collection(artifact2, is_new_version_of=collection)
    assert collection.version == "1"
    assert new_collection.stem_uid == collection.stem_uid
    assert new_collection.version == "2"
    assert new_collection.name == collection.name

    collection.delete(permanent=True)
    artifact2.delete(permanent=True)
    artifact.delete(permanent=True)
