import re

import anndata as ad
import bionty as bt
import lamindb as ln
import numpy as np
import pandas as pd
import pytest
from lamindb.errors import FieldValidationError
from scipy.sparse import csc_matrix, csr_matrix


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
        ln.Collection(artifact, key="Test")
    assert str(error.exconly()).startswith(
        "ValueError: Not all artifacts are yet saved, please save them"
    )
    artifact.save()
    with pytest.raises(ValueError) as error:
        ln.Collection(artifact, artifact)
    assert str(error.exconly()).startswith(
        "ValueError: Only one non-keyword arg allowed: artifacts"
    )
    transform = ln.Transform(key="My test transform").save()
    run = ln.Run(transform).save()
    collection = ln.Collection(artifact, key="My new collection", run=run).save()
    assert collection.run.input_artifacts.get() == artifact
    collection.delete(permanent=True)
    artifact.delete(permanent=True)
    assert ln.Artifact.filter(id=artifact.id).one_or_none() is None


def test_edge_cases(df):
    with pytest.raises(
        FieldValidationError,
        match=re.escape(
            "Only artifacts, key, description, meta, reference, reference_type, run, revises can be passed"
        ),
    ) as error:
        ln.Collection(df, invalid_param=1)

    with pytest.raises(ValueError) as error:
        ln.Collection(1, key="Invalid")
    assert str(error.exconly()).startswith(
        "ValueError: Artifact or list[Artifact] is allowed."
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
    artifact.delete(permanent=True)


def test_from_inconsistent_artifacts(df, adata):
    artifact1 = ln.Artifact.from_df(df, description="My test").save()
    artifact2 = ln.Artifact.from_anndata(adata, description="My test2").save()
    collection = ln.Collection([artifact1, artifact2], name="Inconsistent").save()
    # test idempotency of .save()
    collection.save()
    # create a run context
    ln.track(transform=ln.Transform(key="My test transform"))
    # can iterate over them
    collection.cache()
    assert set(ln.context.run.input_collections.all()) == {collection}
    # loading will throw an error here
    with pytest.raises(ValueError) as error:
        collection.load()
    assert str(error.exconly()).startswith(
        "ValueError: Can only load collections where all artifacts have the same suffix"
    )
    # test through query set
    with pytest.raises(ValueError) as error:
        collection.artifacts.all().load()
    assert str(error.exconly()).startswith(
        "ValueError: Can only load collections where all artifacts have the same suffix"
    )
    collection.describe()
    collection.delete(permanent=True)
    artifact1.delete(permanent=True)
    artifact2.delete(permanent=True)
    ln.context._run = None


def test_from_consistent_artifacts(adata, adata2):
    artifact1 = ln.Artifact.from_anndata(adata, key="my_test.h5ad").save()
    artifact2 = ln.Artifact.from_anndata(adata2, key="my_test.h5ad").save()
    transform = ln.Transform(key="My test transform").save()
    run = ln.Run(transform).save()
    collection = ln.Collection([artifact1, artifact2], key="My test", run=run)
    assert collection._state.adding
    collection.save()
    assert set(collection.run.input_artifacts.all()) == {artifact1, artifact2}
    adata_joined = collection.load()
    assert "artifact_uid" in adata_joined.obs.columns
    assert artifact1.uid in adata_joined.obs.artifact_uid.cat.categories
    # test from query set through collection
    adata_joined = collection.artifacts.order_by("-created_at").load()
    assert "artifact_uid" in adata_joined.obs.columns
    assert artifact1.uid in adata_joined.obs.artifact_uid.cat.categories

    # re-run with hash-based lookup
    collection2 = ln.Collection([artifact1, artifact2], key="My test 1", run=run)
    assert not collection2._state.adding
    assert collection2.id == collection.id
    assert collection2.key == "My test 1"

    collection.delete(permanent=True)
    artifact1.delete(permanent=True)
    artifact2.delete(permanent=True)


def test_mapped(adata, adata2):
    # prepare test data
    adata.strings_to_categoricals()
    adata.obs["feat2"] = adata.obs["feat1"]
    adata.layers["layer1"] = adata.X.copy()
    adata.layers["layer1"][0, 0] = 0
    artifact1 = ln.Artifact.from_anndata(adata, key="part_one.h5ad").save()
    adata2.X = csr_matrix(adata2.X)
    adata2.layers["layer1"] = adata2.X.copy()
    adata2.obs["feat2"] = adata2.obs["feat1"]
    artifact2 = ln.Artifact.from_anndata(
        adata2, key="part_two.zarr", format="zarr"
    ).save()
    adata3 = adata2.copy()
    adata3.var_names = ["A", "B", "C"]
    adata3.obs.loc["0", "feat1"] = np.nan
    artifact3 = ln.Artifact.from_anndata(adata3, key="other_vars.h5ad").save()
    adata4 = adata.copy()
    adata4.layers["layer1"] = csc_matrix(adata4.layers["layer1"])
    artifact4 = ln.Artifact.from_anndata(adata4, description="csc layer").save()
    collection_outer = ln.Collection(
        [artifact1, artifact2, artifact3], key="gather_outer"
    ).save()
    collection_csc = ln.Collection([artifact4, artifact2], key="check_csc").save()
    collection = ln.Collection([artifact1, artifact2], key="gather")
    # test mapped without saving first
    with collection.mapped() as ls_ds:
        assert ls_ds.__class__.__name__ == "MappedCollection"
    collection.save()

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
        # categories in the encoder are sorted
        A_enc = ls_ds.encoders["feat2"]["A"]
        assert A_enc == 0
        B_enc = ls_ds.encoders["feat2"]["B"]
        assert B_enc == 1
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
    assert len(weights) == 4
    assert all(weights == 0.5)
    weights = ls_ds.get_label_weights(["feat1", "feat2"])
    assert len(weights) == 4
    assert all(weights == 0.5)
    weights = ls_ds.get_label_weights(["feat1", "feat2"], scaler=1.0)
    assert all(weights == 1.0 / 3.0)
    weights = ls_ds.get_label_weights(
        ["feat1", "feat2"], scaler=1.0, return_categories=True
    )
    assert weights["A__A"] == 1.0 / 3.0
    assert weights["B__B"] == 1.0 / 3.0

    assert not ls_ds.check_vars_sorted(ascending=True)
    assert not ls_ds.check_vars_sorted(ascending=False)
    assert ls_ds.check_vars_non_aligned(["MYC", "TCF7", "GATA1"]) == []
    ls_ds.var_list = None
    assert not ls_ds.check_vars_sorted()
    ls_ds.var_list = None
    assert ls_ds.check_vars_non_aligned(["MYC", "TCF7", "GATA1"]) == []

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
    # test with QuerySet
    query_set = ln.Artifact.filter(key__in=["part_one.h5ad", "part_two.zarr"])
    with query_set.mapped() as ls_ds:
        assert ls_ds.shape == (4, 3)
    with query_set.order_by("created_at").mapped(stream=True) as ls_ds:
        assert ls_ds.shape == (4, 3)

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
        ls_ds_idx = ls_ds[4]
        assert np.array_equal(ls_ds_idx["X"], np.array([1, 2, 5, 0, 0, 0]))
        assert ls_ds_idx["feat1"] is np.nan
        assert np.array_equal(ls_ds[5]["X"], np.array([4, 5, 8, 0, 0, 0]))
        assert np.issubdtype(ls_ds[2]["X"].dtype, np.integer)
        assert np.issubdtype(ls_ds[4]["X"].dtype, np.integer)
        assert np.array_equal(ls_ds[3]["obsm_X_pca"], np.array([3, 4]))
        assert ls_ds.check_vars_non_aligned(["MYC", "TCF7", "GATA1"]) == [2]
        assert not ls_ds.check_vars_sorted()
        assert len(ls_ds.get_label_weights("feat1")) == 6

    with collection_outer.mapped(layers_keys="layer1", join="outer") as ls_ds:
        assert np.array_equal(ls_ds[0]["layer1"], np.array([0, 0, 0, 3, 0, 2]))
        assert np.array_equal(ls_ds[4]["layer1"], np.array([1, 2, 5, 0, 0, 0]))

    # csc matrix in layers
    with pytest.raises(ValueError):
        collection_csc.mapped(layers_keys="layer1")

    # test with obs_filter
    # tuple as obs_filter is deprecated, test anyways for now
    with collection.mapped(obs_filter=("feat1", ("A", "B"))) as ls_ds:
        assert ls_ds.shape == (4, 3)
        assert np.array_equal(ls_ds[1]["X"], np.array([4, 5, 6]))
        assert np.array_equal(ls_ds[3]["X"], np.array([4, 5, 8]))
        weights = ls_ds.get_label_weights("feat1")
        assert len(weights) == 4
        assert all(weights == 0.5)
    # tuple as obs_filter is deprecated, test anyways for now
    with collection.mapped(obs_filter=("feat1", "B")) as ls_ds:
        assert ls_ds.shape == (2, 3)
        assert np.array_equal(ls_ds[0]["X"], np.array([4, 5, 6]))
        assert np.array_equal(ls_ds[1]["X"], np.array([4, 5, 8]))
        weights = ls_ds.get_label_weights("feat2")
        assert len(weights) == 2
        assert all(weights == 0.5)

    with collection.mapped(obs_filter={"feat1": "B", "feat2": ("A", "B")}) as ls_ds:
        assert ls_ds.shape == (2, 3)
        assert ls_ds.original_shapes == [(1, 3), (1, 3)]
        assert np.array_equal(ls_ds[0]["X"], np.array([4, 5, 6]))
        assert np.array_equal(ls_ds[1]["X"], np.array([4, 5, 8]))
        weights = ls_ds.get_label_weights("feat2")
        assert len(weights) == 2
        assert all(weights == 0.5)
    # nan in filtering values
    with collection_outer.mapped(obs_filter={"feat1": np.nan}, join="outer") as ls_ds:
        assert ls_ds.shape == (1, 6)
        assert np.array_equal(ls_ds[0]["X"], np.array([1, 2, 5, 0, 0, 0]))
    with collection_outer.mapped(
        obs_filter={"feat1": (np.nan,), "feat2": ["A", "B"]}, join="outer"
    ) as ls_ds:
        assert ls_ds.shape == (1, 6)
    with collection_outer.mapped(
        obs_filter={"feat1": (np.nan, "A", "B")}, join="outer"
    ) as ls_ds:
        assert ls_ds.shape == (6, 6)
    with collection_outer.mapped(
        obs_filter={"feat1": ["A", "B"]}, join="outer"
    ) as ls_ds:
        assert ls_ds.shape == (5, 6)
    with collection_outer.mapped(
        obs_filter={"feat1": ("A", np.nan)}, join="outer"
    ) as ls_ds:
        assert ls_ds.shape == (3, 6)

    collection.delete(permanent=True)
    collection_outer.delete(permanent=True)
    collection_csc.delete(permanent=True)
    artifact1.delete(permanent=True)
    artifact2.delete(permanent=True)
    artifact3.delete(permanent=True)
    artifact4.delete(permanent=True)


def test_revise_collection(df, adata):
    # create a versioned collection
    artifact = ln.Artifact.from_df(df, description="test").save()
    collection = ln.Collection(artifact, key="test-collection", version="1")
    assert collection.version == "1"
    assert collection.uid.endswith("0000")
    collection.save()

    artifact = ln.Artifact.from_anndata(adata, description="test").save()

    with pytest.raises(ValueError) as error:
        collection_r2 = ln.Collection(artifact, revises=collection, version="1")
    assert error.exconly() == "ValueError: Please increment the previous version: '1'"

    with pytest.raises(TypeError):
        ln.Collection(adata, revises="wrong-type")

    # create new collection from old collection
    collection_r2 = ln.Collection(artifact, key="test-collection")
    assert collection_r2.stem_uid == collection.stem_uid
    assert collection_r2.uid.endswith("0001")
    # repeat
    collection_r2 = ln.Collection(artifact, key="test-collection")
    assert collection_r2.stem_uid == collection.stem_uid
    assert collection_r2.uid.endswith("0001")
    assert collection_r2.version is None
    assert collection_r2.key == "test-collection"

    collection_r2.save()

    # create new collection from newly versioned collection
    df.iloc[0, 0] = 0
    artifact = ln.Artifact.from_df(df, description="test")
    artifact.save()
    collection_r3 = ln.Collection(
        artifact, key="test-collection", description="test description3", version="2"
    )
    assert collection_r3.stem_uid == collection.stem_uid
    assert collection_r3.version == "2"
    assert collection_r3.uid.endswith("0002")
    assert collection_r3.key == "test-collection"
    assert collection_r3.description == "test description3"

    artifacts_r2 = collection_r2.artifacts.all()
    collection_r2.delete(permanent=True)
    artifacts_r2.delete(permanent=True)
    artifacts = collection.artifacts.all()
    collection.delete(permanent=True)
    artifacts.delete(permanent=True)


def test_collection_append(df, adata):
    artifact = ln.Artifact.from_df(df, description="test").save()
    artifact_1 = ln.Artifact.from_anndata(adata, description="test").save()
    collection = ln.Collection(artifact, key="Test", description="Test append").save()
    new_collection = collection.append(artifact_1).save()

    assert new_collection.key == collection.key
    assert new_collection.description == collection.description
    assert new_collection.uid.endswith("0001")
    artifacts = new_collection.artifacts.all()
    assert len(artifacts) == 2

    new_collection.versions.delete(permanent=True)
    artifacts.delete(permanent=True)


def test_with_metadata(df, adata):
    meta_artifact = ln.Artifact.from_df(df, description="test")
    meta_artifact.save()
    data_artifact = ln.Artifact.from_anndata(adata, description="test adata")
    data_artifact.save()
    collection = ln.Collection(
        data_artifact, name="test collection", meta_artifact=meta_artifact
    )
    collection.save()

    assert collection.meta_artifact == meta_artifact
    assert collection.data_artifact == data_artifact
    collection.delete(permanent=True)
    data_artifact.delete(permanent=True)
    meta_artifact.delete(permanent=True)


def test_collection_get_tracking(df):
    artifact = ln.Artifact.from_df(df, key="df.parquet").save()
    collection = ln.Collection(artifact, key="track-collection").save()

    transform = ln.Transform(key="test track collection via get").save()
    run = ln.Run(transform).save()

    assert (
        ln.Collection.get(key="track-collection", is_run_input=run)
        in run.input_collections.all()
    )

    collection.delete(permanent=True)
    artifact.delete(permanent=True)
    transform.delete()
