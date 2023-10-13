from inspect import signature
from pathlib import Path

import anndata as ad
import lnschema_bionty as lb
import numpy as np
import pandas as pd
import pytest
from django.db.models.deletion import ProtectedError

import lamindb as ln
from lamindb import _dataset

df = pd.DataFrame({"feat1": [1, 2], "feat2": [3, 4]})

adata = ad.AnnData(
    X=np.array([[1, 2, 3], [4, 5, 6]]),
    obs=dict(feat1=["A", "B"]),
    var=pd.DataFrame(index=["MYC", "TCF7", "GATA1"]),
    obsm=dict(X_pca=np.array([[1, 2], [3, 4]])),
)

adata2 = ad.AnnData(
    X=np.array([[1, 2, 5], [4, 5, 8]]),
    obs=dict(feat1=["A", "B"]),
    var=pd.DataFrame(index=["MYC", "TCF7", "GATA1"]),
    obsm=dict(X_pca=np.array([[1, 2], [3, 4]])),
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
        setattr(Mock, name, getattr(_dataset, name))
        assert signature(getattr(Mock, name)) == _dataset.SIGS.pop(name)
    # methods
    for name, sig in _dataset.SIGS.items():
        assert signature(getattr(_dataset, name)) == sig


def test_create_delete_from_single_dataframe():
    df = ln.dev.datasets.df_iris_in_meter_study1()

    dataset = ln.Dataset.from_df(df, name="Iris flower dataset1")
    # because features weren't registered, there is no linked feature set
    assert dataset._feature_sets == {}

    # register features
    ln.save(ln.Feature.from_df(df))

    # won't work with features like so
    dataset = ln.Dataset(df, name="Iris flower dataset1")
    assert dataset._feature_sets == {}

    # will work like so
    dataset = ln.Dataset.from_df(df, name="Iris flower dataset1")
    assert "columns" in dataset._feature_sets

    dataset.save()

    # basics
    assert dataset.load().iloc[0].tolist() == df.iloc[0].tolist()
    file = dataset.file
    assert file.description == f"See dataset {dataset.uid}"
    assert dataset.hash == file.hash
    assert dataset.uid == file.uid
    assert ln.File.filter(uid=dataset.uid).one_or_none() is not None
    assert ln.File.filter(uid=file.uid).one_or_none() is not None

    # features
    feature_list = [
        "sepal_length",
        "sepal_width",
        "petal_length",
        "petal_width",
        "iris_species_name",
    ]
    assert len(ln.Feature.filter(name__in=feature_list).list()) == 5
    feature_set = ln.FeatureSet.filter(datasets=dataset).one()
    feature_list_queried = ln.Feature.filter(feature_sets=feature_set).list()
    feature_list_queried = [feature.name for feature in feature_list_queried]
    assert set(feature_list_queried) == set(feature_list)
    # the feature_set is also linked to the file
    assert ln.FeatureSet.filter(files=dataset.file).one() == feature_set

    # accidental recreation (re-load based on hash)
    dataset1 = ln.Dataset.from_df(df, name="Iris Flower data1")
    assert dataset1.id == dataset.id
    assert dataset1.hash == dataset.hash

    # now proceed to deletion
    dataset.delete(storage=True)
    assert ln.File.filter(uid=dataset.uid).one_or_none() is None
    assert ln.File.filter(uid=file.uid).one_or_none() is None


def test_create_delete_from_single_anndata():
    ln.track(ln.Transform(name="Test transform"))
    dataset = ln.Dataset(adata, name="My adata")
    dataset.save()
    dataset.delete(storage=True)
    assert ln.File.filter(id=dataset.id).one_or_none() is None
    assert ln.File.filter(id=dataset.file.id).one_or_none() is None
    # and now with from_anndata
    lb.settings.species = "human"
    dataset = ln.Dataset.from_anndata(adata, name="My adata", field=lb.Gene.symbol)
    # let's now try passing an AnnData-like file with some feature sets linked
    ln.save(ln.Feature.from_df(adata.obs))
    file = ln.File.from_anndata(adata, description="my adata", field=lb.Gene.symbol)
    file.save()
    ln.save(lb.Gene.from_values(adata.var.index, "symbol"))
    dataset = ln.Dataset.from_anndata(file, name="My dataset", field=lb.Gene.symbol)
    dataset.save()
    dataset.describe()
    dataset.view_flow()
    feature_sets_queried = dataset.feature_sets.all()
    features_queried = ln.Feature.filter(feature_sets__in=feature_sets_queried).all()
    assert set(features_queried.list("name")) == set(adata.obs.columns)
    genes_queried = lb.Gene.filter(feature_sets__in=feature_sets_queried).all()
    assert set(genes_queried.list("symbol")) == set(adata.var.index)
    feature_sets_queried.delete()
    features_queried.delete()
    genes_queried.delete()
    dataset.delete(storage=True)
    ln.dev.run_context.run = None
    ln.dev.run_context.transform = None


def test_from_single_file():
    lb.settings.species = "human"
    features = ln.Feature.from_df(adata.obs)
    validated = ln.Feature.validate(
        [feature.name for feature in features], field="name"
    )
    ln.save([feature for (feature, valid) in zip(features, validated) if valid])
    file = ln.File.from_anndata(adata, description="My adata", field=lb.Gene.symbol)
    with pytest.raises(ValueError) as error:
        ln.Dataset(file)
    assert str(error.exconly()).startswith(
        "ValueError: Save file before creating dataset!"
    )
    file.save()
    with pytest.raises(ValueError) as error:
        ln.Dataset(file, file)
    assert str(error.exconly()).startswith(
        "ValueError: Only one non-keyword arg allowed: data"
    )
    transform = ln.Transform(name="My test transform")
    transform.save()
    run = ln.Run(transform)
    run.save()
    dataset = ln.Dataset(file, name="My new dataset", run=run)
    dataset.save()
    # test data flow
    assert dataset.run.input_files.get() == file
    # test features
    assert set(file.feature_sets.list("id")) == set(
        dataset.file.feature_sets.list("id")
    )
    assert set(file.features._feature_set_by_slot.keys()) == set(
        dataset.features._feature_set_by_slot.keys()
    )
    feature_sets_queried = file.feature_sets.all()
    features_queried = ln.Feature.filter(feature_sets__in=feature_sets_queried).all()
    feature_sets_queried.delete()
    features_queried.delete()
    dataset.delete(storage=True)
    assert ln.File.filter(id=dataset.id).one_or_none() is None
    assert ln.File.filter(id=dataset.file.id).one_or_none() is None


def test_edge_cases():
    with pytest.raises(ValueError) as error:
        ln.Dataset(df, invalid_param=1)
    assert str(error.exconly()).startswith(
        "ValueError: Only data, name, run, description, reference, reference_type can be passed, you passed: "  # noqa
    )
    with pytest.raises(ValueError) as error:
        ln.Dataset(1, name="Invalid")
    assert str(error.exconly()).startswith(
        "ValueError: Only DataFrame, AnnData and iterable of File is allowed"
    )
    file = ln.File(df, description="Test file")
    assert file._state.adding
    with pytest.raises(ValueError) as error:
        ln.Dataset([file])
    assert str(error.exconly()).startswith(
        "ValueError: Not all files are yet saved, please save them"
    )
    file.save()
    with pytest.raises(ValueError) as error:
        ln.Dataset([file, file])
    assert str(error.exconly()).startswith(
        "ValueError: Please pass files with distinct hashes: these ones are non-unique"
    )
    file.delete(storage=True)


def test_backed():
    dataset = ln.Dataset(adata, name="My test")
    dataset.backed()


def test_from_inconsistent_files():
    file1 = ln.File(df, description="My test")
    file1.save()
    file2 = ln.File(adata, description="My test2")
    file2.save()
    dataset = ln.Dataset([file1, file2], name="Inconsistent")
    dataset.save()
    # create a run context
    ln.track(ln.Transform(name="My test transform"))
    # can iterate over them
    files = dataset.files.all()  # noqa
    assert set(ln.dev.run_context.run.input_datasets.all()) == {dataset}
    # loading will throw an error here
    with pytest.raises(RuntimeError) as error:
        dataset.load()
    assert str(error.exconly()).startswith(
        "RuntimeError: Can only load datasets where all files have the same suffix"
    )
    file1.delete(storage=True)
    file2.delete(storage=True)
    dataset.delete()
    ln.dev.run_context.run = None
    ln.dev.run_context.transform = None


def test_from_consistent_files():
    file1 = ln.File(adata, description="My test")
    file1.save()
    file2 = ln.File(adata2, description="My test2")
    file2.save()
    transform = ln.Transform(name="My test transform")
    transform.save()
    run = ln.Run(transform)
    run.save()
    dataset = ln.Dataset([file1, file2], name="My test", run=run)
    dataset.save()
    assert set(dataset.run.input_files.all()) == {file1, file2}
    dataset.load()
    with pytest.raises(RuntimeError) as error:
        dataset.backed()
    assert str(error.exconly()).startswith(
        "RuntimeError: Can only call backed() for datasets with a single file"
    )
    file1.delete(storage=True)
    file2.delete(storage=True)
    dataset.delete()


def test_is_new_version_of_versioned_dataset():
    # create a versioned dataset
    dataset = ln.Dataset(df, name="test", version="1")
    assert dataset.version == "1"
    dataset.save()

    with pytest.raises(ValueError) as error:
        dataset_v2 = ln.Dataset(adata, is_new_version_of=dataset, version="1")
    assert error.exconly() == "ValueError: Please increment the previous version: '1'"

    # create new dataset from old dataset
    dataset_v2 = ln.Dataset(adata, is_new_version_of=dataset)
    assert dataset.version == "1"
    assert (
        dataset.initial_version_id is None
    )  # initial dataset has initial_version_id None
    assert dataset_v2.initial_version_id == dataset.id
    assert dataset_v2.version == "2"
    assert dataset_v2.name == "test"

    dataset_v2.save()

    # create new dataset from newly versioned dataset
    df.iloc[0, 0] = 0
    dataset_v3 = ln.Dataset(df, name="test1", is_new_version_of=dataset_v2)
    assert dataset_v3.initial_version_id == dataset.id
    assert dataset_v3.version == "3"
    assert dataset_v3.name == "test1"

    # test that reference dataset cannot be deleted
    with pytest.raises(ProtectedError):
        dataset.delete(storage=True)
    dataset_v2.delete(storage=True)
    dataset.delete(storage=True)


def test_is_new_version_of_unversioned_dataset():
    # unversioned dataset
    dataset = ln.Dataset(df, name="test2")
    assert dataset.initial_version_id is None
    assert dataset.version is None

    # what happens if we don't save the old dataset?
    # add a test for it!
    dataset.save()

    with pytest.raises(TypeError):
        ln.Dataset(adata, is_new_version_of="wrong-type")

    # create new dataset from old dataset
    new_dataset = ln.Dataset(adata, is_new_version_of=dataset)
    assert dataset.version == "1"
    assert dataset.initial_version is None
    assert new_dataset.initial_version_id == dataset.id
    assert new_dataset.version == "2"
    assert new_dataset.name == dataset.name

    dataset.delete(storage=True)


def test_dataset_from_storage():
    # this should only accept folders
    path = Path("./random_storage/")
    path.mkdir(exist_ok=True)
    dataset = ln.Dataset(path, name="My test dataset")
    assert dataset.path.name == "random_storage"
    dataset.save()
    dataset.delete()
