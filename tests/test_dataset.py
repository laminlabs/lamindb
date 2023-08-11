import anndata as ad
import numpy as np

import lamindb as ln

adata = ad.AnnData(
    X=np.array([[1, 2, 3], [4, 5, 6]]),
    obs=dict(Obs=["A", "B"]),
    var=dict(Feat=["a", "b", "c"]),
    obsm=dict(X_pca=np.array([[1, 2], [3, 4]])),
)


def test_create_delete_from_single_dataframe():
    df = ln.dev.datasets.df_iris_in_meter_batch1()

    # try to create and save the dataset
    dataset = ln.Dataset(df, name="Iris flower dataset1")
    # because features weren't registered, there is no linked feature set
    assert dataset._feature_sets == {}
    # register features and then repeat
    ln.save(ln.Feature.from_df(df))
    dataset = ln.Dataset(df, name="Iris flower dataset1")
    assert "columns" in dataset._feature_sets

    dataset.save()

    # basics
    assert dataset.load().iloc[0].tolist() == df.iloc[0].tolist()
    file = dataset.file
    assert file.description == f"See dataset {dataset.id}"
    assert dataset.hash == file.hash
    assert dataset.id == file.id
    assert ln.File.filter(id=dataset.id).one_or_none() is not None
    assert ln.File.filter(id=file.id).one_or_none() is not None

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

    # now proceed to deletion
    dataset.delete(storage=True)
    assert ln.File.filter(id=dataset.id).one_or_none() is None
    assert ln.File.filter(id=file.id).one_or_none() is None


# def test_create_delete_from_single_anndata():
#     dataset = ln.Dataset(adata, name="My adata")
#     dataset.save()
#     assert dataset.load().iloc[0].tolist() == df.iloc[0].tolist()
#     file = dataset.file
#     assert file.description is None
#     assert dataset.hash == file.hash
#     assert dataset.id == file.id
#     assert ln.File.filter(id=dataset.id).one_or_none() is not None
#     assert ln.File.filter(id=file.id).one_or_none() is not None
#     dataset.delete(storage=True)
#     assert ln.File.filter(id=dataset.id).one_or_none() is None
#     assert ln.File.filter(id=file.id).one_or_none() is None
