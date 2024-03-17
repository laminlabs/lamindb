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
    class_methods = ["from_df", "from_anndata"]
    for name in class_methods:
        setattr(Mock, name, getattr(_collection, name))
        assert signature(getattr(Mock, name)) == _collection.SIGS.pop(name)
    # methods
    for name, sig in _collection.SIGS.items():
        assert signature(getattr(_collection, name)) == sig


def test_create_delete_from_single_dataframe(df):
    df = ln.core.datasets.df_iris_in_meter_study1()

    collection = ln.Collection.from_df(df, name="Iris flower collection1")
    collection.save()

    # register features
    ln.save(ln.Feature.from_df(df))

    # link features to collection
    features = ln.Feature.from_df(df)
    collection.features.add(features)
    assert collection.features["columns"] is not None

    # basics
    assert collection.load().iloc[0].tolist() == df.iloc[0].tolist()
    artifact = collection.artifact
    assert artifact.description == f"See collection {collection.uid}"
    assert collection.hash == artifact.hash
    assert collection.uid == artifact.uid
    assert ln.Artifact.filter(uid=collection.uid).one_or_none() is not None
    assert ln.Artifact.filter(uid=artifact.uid).one_or_none() is not None

    # features
    feature_list = [
        "sepal_length",
        "sepal_width",
        "petal_length",
        "petal_width",
        "iris_organism_name",
    ]
    assert len(ln.Feature.filter(name__in=feature_list).list()) == 5
    feature_set = ln.FeatureSet.filter(collections=collection).one()
    feature_list_queried = ln.Feature.filter(feature_sets=feature_set).list()
    feature_list_queried = [feature.name for feature in feature_list_queried]
    assert set(feature_list_queried) == set(feature_list)
    # # the feature_set is also linked to the artifact
    # assert ln.FeatureSet.filter(artifacts=collection.artifact).one() == feature_set

    # accidental recreation (re-load based on hash)
    collection1 = ln.Collection.from_df(df, name="Iris Flower data1")
    assert collection1.id == collection.id
    assert collection1.hash == collection.hash

    # now proceed to deletion
    collection.delete(permanent=True)
    collection.artifact.delete(permanent=True, storage=True)
    assert ln.Artifact.filter(uid=collection.uid).one_or_none() is None
    assert ln.Artifact.filter(uid=artifact.uid).one_or_none() is None


def test_create_delete_from_single_anndata(adata):
    ln.track(transform=ln.Transform(name="Test transform"))
    collection = ln.Collection.from_anndata(adata, name="My adata")
    collection.save()
    collection.delete(permanent=True)
    collection.artifact.delete(permanent=True, storage=True)
    assert ln.Artifact.filter(id=collection.id).one_or_none() is None
    assert ln.Artifact.filter(id=collection.artifact.id).one_or_none() is None
    # and now with from_anndata
    bt.settings.organism = "human"
    collection = ln.Collection.from_anndata(adata, name="My adata")
    # let's now try passing an AnnData-like artifact with some feature sets linked
    ln.save(ln.Feature.from_df(adata.obs))
    artifact = ln.Artifact.from_anndata(adata, description="my adata")
    artifact.save()
    ln.save(bt.Gene.from_values(adata.var.index, "symbol"))
    collection = ln.Collection.from_anndata(artifact, name="My collection")
    collection.save()
    collection.describe()
    collection.view_lineage()
    collection.features.add_from_anndata(var_field=bt.Gene.symbol)
    feature_sets_queried = collection.feature_sets.all()
    features_queried = ln.Feature.filter(feature_sets__in=feature_sets_queried).all()
    assert set(features_queried.list("name")) == set(adata.obs.columns)
    genes_queried = bt.Gene.filter(feature_sets__in=feature_sets_queried).all()
    assert set(genes_queried.list("symbol")) == set(adata.var.index)
    feature_sets_queried.delete()
    features_queried.delete()
    genes_queried.delete()
    collection.delete(permanent=True)
    collection.artifact.delete(permanent=True, storage=True)
    ln.core.run_context.run = None
    ln.core.run_context.transform = None


def test_from_single_artifact(adata):
    bt.settings.organism = "human"
    features = ln.Feature.from_df(adata.obs)
    validated = ln.Feature.validate(
        [feature.name for feature in features], field="name"
    )
    ln.save([feature for (feature, valid) in zip(features, validated) if valid])
    artifact = ln.Artifact.from_anndata(adata, description="My adata")
    with pytest.raises(ValueError) as error:
        ln.Collection(artifact)
    assert str(error.exconly()).startswith(
        "ValueError: Save artifact before creating collection!"
    )
    artifact.save()
    with pytest.raises(ValueError) as error:
        ln.Collection(artifact, artifact)
    assert str(error.exconly()).startswith(
        "ValueError: Only one non-keyword arg allowed: data"
    )
    transform = ln.Transform(name="My test transform")
    transform.save()
    run = ln.Run(transform)
    run.save()
    collection = ln.Collection(artifact, name="My new collection", run=run)
    collection.save()
    # test data flow
    assert collection.run.input_artifacts.get() == artifact
    # test features
    artifact.features.add_from_anndata(var_field=bt.Gene.symbol)
    collection.features.add_from_anndata(var_field=bt.Gene.symbol)
    assert set(artifact.feature_sets.list("id")) == set(
        collection.artifact.feature_sets.list("id")
    )
    assert set(artifact.features._feature_set_by_slot.keys()) == set(
        collection.features._feature_set_by_slot.keys()
    )
    feature_sets_queried = artifact.feature_sets.all()
    features_queried = ln.Feature.filter(feature_sets__in=feature_sets_queried).all()
    feature_sets_queried.delete()
    features_queried.delete()
    collection.delete(permanent=True)
    collection.artifact.delete(permanent=True, storage=True)
    assert ln.Artifact.filter(id=collection.id).one_or_none() is None
    assert ln.Artifact.filter(id=collection.artifact.id).one_or_none() is None


def test_edge_cases(df):
    with pytest.raises(ValueError) as error:
        ln.Collection(df, invalid_param=1)
    assert str(error.exconly()).startswith(
        "ValueError: Only data, name, run, description, reference, reference_type, visibility can be passed, you passed: "
    )
    with pytest.raises(ValueError) as error:
        ln.Collection(1, name="Invalid")
    assert str(error.exconly()).startswith(
        "ValueError: Only DataFrame, AnnData, Artifact or list of artifacts is allowed."
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


def test_backed(adata):
    collection = ln.Collection.from_anndata(adata, name="My test")
    collection.backed()


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
    artifacts = collection.artifacts.all()  # noqa
    assert set(ln.core.run_context.run.input_collections.all()) == {collection}
    # loading will throw an error here
    with pytest.raises(RuntimeError) as error:
        collection.load()
    assert str(error.exconly()).startswith(
        "RuntimeError: Can only load collections where all artifacts have the same suffix"
    )
    artifact1.delete(permanent=True, storage=True)
    artifact2.delete(permanent=True, storage=True)
    collection.delete(permanent=True)
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
    collection.save()
    assert set(collection.run.input_artifacts.all()) == {artifact1, artifact2}
    adata_joined = collection.load()
    assert "artifact_uid" in adata_joined.obs.columns
    assert artifact1.uid in adata_joined.obs.artifact_uid.cat.categories
    with pytest.raises(RuntimeError) as error:
        collection.backed()
    assert str(error.exconly()).startswith(
        "RuntimeError: Can only call backed() for collections with a single artifact"
    )
    artifact1.delete(permanent=True, storage=True)
    artifact2.delete(permanent=True, storage=True)
    collection.delete(permanent=True)


def test_collection_mapped(adata, adata2):
    adata.strings_to_categoricals()
    adata.obs["feat2"] = adata.obs["feat1"]
    artifact1 = ln.Artifact.from_anndata(adata, description="Part one")
    artifact1.save()
    adata2.X = csr_matrix(adata2.X)
    adata2.obs["feat2"] = adata2.obs["feat1"]
    artifact2 = ln.Artifact.from_anndata(adata2, description="Part two", format="zrad")
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
        ls_ds = collection.mapped(label_keys="feat1", encode_labels=["feat3"])
    with pytest.raises(ValueError):
        ls_ds = collection.mapped(
            label_keys="feat1", unknown_label={"feat3": "Unknown"}
        )
    with collection.mapped(label_keys=["feat1", "feat2"], unknown_label="A") as ls_ds:
        assert ls_ds.encoders["feat1"]["A"] == -1
        assert ls_ds.encoders["feat1"]["B"] == 0
        assert ls_ds.encoders["feat2"]["A"] == -1
        assert ls_ds.encoders["feat2"]["B"] == 0
        assert ls_ds[0]["feat1"] == -1
        assert ls_ds[1]["feat1"] == 0
        assert ls_ds[0]["feat2"] == -1
        assert ls_ds[1]["feat2"] == 0
    with collection.mapped(
        label_keys=["feat1", "feat2"], unknown_label={"feat1": "A"}
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
        label_keys=["feat1", "feat2"], unknown_label="A", encode_labels=["feat1"]
    ) as ls_ds:
        assert ls_ds.encoders["feat1"]["A"] == -1
        assert ls_ds.encoders["feat1"]["B"] == 0
        assert "feat2" not in ls_ds.encoders
        assert ls_ds[0]["feat1"] == -1
        assert ls_ds[1]["feat1"] == 0
        assert ls_ds[0]["feat2"] == "A"
        assert ls_ds[1]["feat2"] == "B"

    ls_ds = collection.mapped(label_keys="feat1")
    assert not ls_ds.closed

    assert len(ls_ds) == 4
    assert len(ls_ds[0]) == 3 and len(ls_ds[2]) == 3
    assert len(ls_ds[0]["x"]) == 3
    assert np.array_equal(ls_ds[2]["x"], np.array([1, 2, 5]))
    weights = ls_ds.get_label_weights("feat1")
    assert all(weights[1:] == weights[0])
    weights = ls_ds.get_label_weights(["feat1", "feat2"])
    assert all(weights[1:] == weights[0])
    ls_ds.close()
    assert ls_ds.closed
    del ls_ds

    with collection.mapped(label_keys="feat1", join="inner", dtype="float32") as ls_ds:
        assert not ls_ds.closed
        assert len(ls_ds) == 4
        assert len(ls_ds[0]) == 3 and len(ls_ds[2]) == 3
        assert str(ls_ds[0]["x"].dtype) == "float32"
        assert str(ls_ds[2]["x"].dtype) == "float32"
    assert ls_ds.closed

    ls_ds = collection.mapped(label_keys="feat1", parallel=True)
    assert len(ls_ds[0]) == 3 and len(ls_ds[2]) == 3
    assert ls_ds[0]["_storage_idx"] == 0
    assert ls_ds[2]["_storage_idx"] == 1

    with pytest.raises(ValueError):
        with collection_outer.mapped(label_keys="feat1", join="inner"):
            pass

    with collection_outer.mapped(label_keys="feat1", join="outer") as ls_ds:
        assert ls_ds.join_vars == "outer"
        assert len(ls_ds.var_joint) == 6
        assert len(ls_ds[0]) == 3
        assert len(ls_ds[0]["x"]) == 6
        assert np.array_equal(ls_ds[0]["x"], np.array([0, 0, 0, 3, 1, 2]))
        assert np.array_equal(ls_ds[1]["x"], np.array([0, 0, 0, 6, 4, 5]))
        assert np.array_equal(ls_ds[2]["x"], np.array([0, 0, 0, 5, 1, 2]))
        assert np.array_equal(ls_ds[3]["x"], np.array([0, 0, 0, 8, 4, 5]))
        assert np.array_equal(ls_ds[4]["x"], np.array([1, 2, 5, 0, 0, 0]))
        assert np.array_equal(ls_ds[5]["x"], np.array([4, 5, 8, 0, 0, 0]))
        assert np.issubdtype(ls_ds[2]["x"].dtype, np.integer)
        assert np.issubdtype(ls_ds[4]["x"].dtype, np.integer)

    artifact1.delete(permanent=True, storage=True)
    artifact2.delete(permanent=True, storage=True)
    artifact3.delete(permanent=True, storage=True)
    collection.delete(permanent=True)
    collection_outer.delete(permanent=True)


def test_is_new_version_of_versioned_collection(df, adata):
    # create a versioned collection
    collection = ln.Collection.from_df(df, name="test", version="1")
    assert collection.version == "1"
    collection.save()

    with pytest.raises(ValueError) as error:
        collection_v2 = ln.Collection.from_anndata(
            adata, is_new_version_of=collection, version="1"
        )
    assert error.exconly() == "ValueError: Please increment the previous version: '1'"

    # create new collection from old collection
    collection_v2 = ln.Collection.from_anndata(adata, is_new_version_of=collection)
    assert collection.version == "1"
    assert collection_v2.stem_uid == collection.stem_uid
    assert collection_v2.version == "2"
    assert collection_v2.name == "test"

    collection_v2.save()

    # create new collection from newly versioned collection
    df.iloc[0, 0] = 0
    collection_v3 = ln.Collection.from_df(
        df, name="test1", is_new_version_of=collection_v2
    )
    assert collection_v3.stem_uid == collection.stem_uid
    assert collection_v3.version == "3"
    assert collection_v3.name == "test1"

    # test that reference collection cannot be deleted
    collection_v2.delete(permanent=True)
    collection_v2.artifact.delete(permanent=True, storage=True)
    collection.delete(permanent=True)
    collection.artifact.delete(permanent=True, storage=True)


def test_is_new_version_of_unversioned_collection(df, adata):
    # unversioned collection
    collection = ln.Collection.from_df(df, name="test2")
    assert collection.version is None

    # what happens if we don't save the old collection?
    # add a test for it!
    collection.save()

    with pytest.raises(TypeError):
        ln.Collection.from_anndata(adata, is_new_version_of="wrong-type")

    # create new collection from old collection
    new_collection = ln.Collection.from_anndata(adata, is_new_version_of=collection)
    assert collection.version == "1"
    assert new_collection.stem_uid == collection.stem_uid
    assert new_collection.version == "2"
    assert new_collection.name == collection.name

    collection.delete(permanent=True)
    collection.artifact.delete(permanent=True, storage=True)
