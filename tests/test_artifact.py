"""Artifact tests.

Also see `test_artifact_folders.py` for tests of folder-like artifacts.

"""

import shutil
from inspect import signature
from pathlib import Path

import anndata as ad
import bionty as bt
import lamindb as ln
import lamindb_setup
import numpy as np
import pandas as pd
import pytest
from lamindb import _artifact
from lamindb._artifact import (
    check_path_is_child_of_root,
    data_is_anndata,
    get_relative_path_to_directory,
    process_data,
)
from lamindb.core._settings import settings
from lamindb.core.exceptions import IntegrityError
from lamindb.core.storage._zarr import write_adata_zarr, zarr_is_adata
from lamindb.core.storage.paths import (
    AUTO_KEY_PREFIX,
    auto_storage_key_from_artifact_uid,
    delete_storage,
    load_to_memory,
    read_fcs,
    read_tsv,
)
from lamindb_setup.core.upath import (
    CloudPath,
    LocalPathClasses,
    UPath,
)

# how do we properly abstract out the default storage variable?
# currently, we're only mocking it through `default_storage` as
# set in conftest.py

ln.settings.verbosity = "success"
bt.settings.organism = "human"


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
def adata_file():
    adata = ad.AnnData(
        X=np.array([[1, 2, 3], [4, 5, 6]]),
        obs={"feat1": ["A", "B"]},
        var=pd.DataFrame(index=["MYC", "TCF7", "GATA1"]),
        obsm={"X_pca": np.array([[1, 2], [3, 4]])},
    )
    filepath = Path("adata_file.h5ad")
    adata.write(filepath)
    yield "adata_file.h5ad"
    filepath.unlink()


@pytest.fixture
def data(request):
    if request.param == "adata":
        return request.getfixturevalue("adata")
    else:
        return request.param


@pytest.fixture(scope="module")
def tsv_file():
    filepath = Path("test.tsv")
    pd.DataFrame([1, 2]).to_csv(filepath, sep="\t")
    yield filepath
    filepath.unlink()


@pytest.fixture(scope="module")
def zip_file():
    filepath = Path("test.zip")
    pd.DataFrame([1, 2]).to_csv(filepath, sep="\t")
    yield filepath
    filepath.unlink()


@pytest.fixture(scope="module")
def fcs_file():
    return ln.core.datasets.file_fcs()


def test_signatures():
    # this seems currently the easiest and most transparent
    # way to test violations of the signature equality
    # the MockORM class is needed to get inspect.signature
    # to work
    class Mock:
        pass

    # class methods
    class_methods = ["from_dir", "from_df", "from_anndata", "from_mudata"]
    for name in class_methods:
        setattr(Mock, name, getattr(_artifact, name))
        assert signature(getattr(Mock, name)) == _artifact.SIGS.pop(name)
    # methods
    for name, sig in _artifact.SIGS.items():
        if name == "delete":
            # have a temporary fix in delete regarding "using_key"
            continue
        assert signature(getattr(_artifact, name)) == sig


def test_data_is_anndata_paths():
    assert data_is_anndata("something.h5ad")
    assert data_is_anndata("something.anndata.zarr")
    assert data_is_anndata("s3://somewhere/something.anndata.zarr")
    assert not data_is_anndata("s3://somewhere/something.zarr")


def test_basic_validation():
    # extra kwargs
    with pytest.raises(ValueError):
        ln.Artifact("testpath.csv", description="test1b", extra_kwarg="extra")

    # > 1 args
    with pytest.raises(ValueError) as error:
        ln.Artifact("testpath.csv", "testpath.csv")
    assert error.exconly() == "ValueError: Only one non-keyword arg allowed: data"


# comment out for now
# AUTO_KEY_PREFIX
#    with pytest.raises(ValueError) as error:
#        ln.Artifact.from_df(df, key=".lamindb/test_df.parquet")
#    assert error.exconly() == "ValueError: Key cannot start with .lamindb/"


def test_revise_artifact(df, adata):
    # attempt to create a file with an invalid version
    with pytest.raises(ValueError) as error:
        artifact = ln.Artifact.from_df(df, description="test", version=0)
    assert (
        error.exconly()
        == "ValueError: `version` parameter must be `None` or `str`, e.g., '0.1', '1',"
        " '2', etc."
    )

    # create a file and tag it with a version
    artifact = ln.Artifact.from_df(df, description="test", version="1")
    assert artifact.version == "1"
    assert artifact.uid.endswith("0000")

    assert artifact.path.exists()  # because of cache file already exists
    artifact.save()
    assert artifact.path.exists()

    with pytest.raises(ValueError) as error:
        artifact_r2 = ln.Artifact.from_anndata(adata, revises=artifact, version="1")
    assert error.exconly() == "ValueError: Please increment the previous version: '1'"

    # create new file from old file
    artifact_r2 = ln.Artifact.from_anndata(adata, revises=artifact)
    assert artifact_r2.uid.endswith("0001")
    assert artifact_r2.stem_uid == artifact.stem_uid
    assert artifact_r2.version is None
    assert artifact_r2.key is None
    assert artifact_r2.description == "test"
    assert artifact_r2._revises is not None
    artifact_r2.save()
    assert artifact_r2.path.exists()
    assert artifact_r2._revises is None

    # create new file from newly versioned file
    df.iloc[0, 0] = 0  # mutate dataframe so that hash lookup doesn't trigger
    artifact_r3 = ln.Artifact.from_df(
        df, description="test1", revises=artifact_r2, version="2"
    )
    assert artifact_r3.uid.endswith("0002")
    assert artifact_r3.stem_uid == artifact.stem_uid
    assert artifact_r3.version == "2"
    assert artifact_r3.description == "test1"

    # revise by matching on `key`
    key = "my-test-dataset.parquet"
    artifact_r2.key = key
    artifact_r2.save()
    artifact_r3 = ln.Artifact.from_df(df, description="test1", key=key, version="2")
    assert artifact_r3.uid.endswith("0002")
    assert artifact_r3.stem_uid == artifact.stem_uid
    assert artifact_r3.key == key
    assert artifact_r3.version == "2"
    assert artifact_r3.description == "test1"

    artifact_r3 = ln.Artifact.from_df(
        df, description="test1", key="my-test-dataset1.parquet", version="2"
    )

    with pytest.raises(TypeError) as error:
        ln.Artifact.from_df(df, description="test1a", revises=ln.Transform())
    assert error.exconly() == "TypeError: `revises` has to be of type `Artifact`"

    artifact_r2.delete(permanent=True, storage=True)
    artifact.delete(permanent=True, storage=True)

    # unversioned file
    artifact = ln.Artifact.from_df(df, description="test2")
    assert artifact.version is None

    # what happens if we don't save the old file?
    # add a test for it!
    artifact.save()

    # create new file from old file
    new_artifact = ln.Artifact.from_anndata(adata, revises=artifact)
    assert artifact.version is None
    assert new_artifact.stem_uid == artifact.stem_uid
    assert new_artifact.version is None
    assert new_artifact.description == artifact.description

    artifact.delete(permanent=True, storage=True)


# also test legacy name parameter (got removed by description)
def test_create_from_dataframe(df):
    artifact = ln.Artifact.from_df(df, description="test1")
    assert artifact.description == "test1"
    assert artifact.key is None
    assert artifact._accessor == "DataFrame"
    assert artifact.type == "dataset"
    assert hasattr(artifact, "_local_filepath")
    artifact.save()
    # can't do backed
    with pytest.raises(ValueError):
        artifact.open()
    # check that the local filepath has been cleared
    assert not hasattr(artifact, "_local_filepath")
    artifact.delete(permanent=True, storage=True)


def test_create_from_dataframe_using_from_df_and_link_features(df):
    description = "my description"
    artifact = ln.Artifact.from_df(
        df, key="folder/hello.parquet", description=description
    )
    with pytest.raises(ValueError):
        artifact.features["columns"]
    artifact = ln.Artifact.from_df(df, description=description)
    # backward compatibility for ln.Artifact to take a DataFrame
    artifact = ln.Artifact(df, key="folder/hello.parquet", description=description)
    assert artifact.description == description
    assert artifact._accessor == "DataFrame"
    assert artifact.key == "folder/hello.parquet"
    assert artifact._key_is_virtual
    assert artifact.uid in artifact.path.as_posix()
    artifact.save()
    # register features from df columns
    features = ln.Feature.from_df(df)
    ln.save(features)
    # link features
    artifact.features._add_set_from_df()
    # mere access test right now
    artifact.features["columns"]
    feature_set_queried = artifact.feature_sets.get()  # exactly one
    feature_list_queried = ln.Feature.filter(feature_sets=feature_set_queried).list()
    feature_list_queried = [feature.name for feature in feature_list_queried]
    assert set(feature_list_queried) == set(df.columns)
    artifact.delete(permanent=True, storage=True)
    feature_set_queried.delete()
    ln.Feature.filter(name__in=["feat1", "feat2"]).delete()


def test_create_from_anndata_in_memory_and_link_features(adata):
    ln.save(
        bt.Gene.from_values(adata.var.index, field=bt.Gene.symbol, organism="human")
    )
    ln.save(ln.Feature.from_df(adata.obs))
    artifact = ln.Artifact.from_anndata(adata, description="test")
    assert artifact._accessor == "AnnData"
    assert hasattr(artifact, "_local_filepath")
    artifact.save()
    # check that the local filepath has been cleared
    assert not hasattr(artifact, "_local_filepath")
    # link features
    artifact.features._add_set_from_anndata(var_field=bt.Gene.symbol, organism="human")
    feature_sets_queried = artifact.feature_sets.all()
    features_queried = ln.Feature.filter(feature_sets__in=feature_sets_queried).all()
    assert set(features_queried.list("name")) == set(adata.obs.columns)
    genes_queried = bt.Gene.filter(feature_sets__in=feature_sets_queried).all()
    assert set(genes_queried.list("symbol")) == set(adata.var.index)
    artifact.delete(permanent=True, storage=True)
    feature_sets_queried.delete()
    features_queried.delete()
    genes_queried.delete()


def test_create_from_anndata_strpath(adata_file):
    artifact = ln.Artifact.from_anndata(adata_file, description="test adata file")
    artifact.save()
    assert artifact._accessor == "AnnData"
    artifact.delete(permanent=True, storage=True)


@pytest.mark.parametrize(
    "data", ["adata", "s3://lamindb-test/scrnaseq_pbmc68k_tiny.h5ad"], indirect=True
)
def test_create_from_anndata_in_storage(data):
    if isinstance(data, ad.AnnData):
        artifact = ln.Artifact.from_anndata(
            data, description="test_create_from_anndata"
        )
        assert artifact._accessor == "AnnData"
        assert hasattr(artifact, "_local_filepath")
    else:
        previous_storage = ln.setup.settings.storage.root_as_str
        ln.settings.storage = "s3://lamindb-test"
        filepath = data
        # TODO: automatically add accessor based on file suffix
        artifact = ln.Artifact(filepath)
    artifact.save()
    # check that the local filepath has been cleared
    assert not hasattr(artifact, "_local_filepath")
    if not isinstance(data, ad.AnnData):
        ln.settings.storage = previous_storage


# this tests the basic (non-provenance-related) metadata
@pytest.mark.parametrize("key_is_virtual", [True, False])
@pytest.mark.parametrize("key", [None, "my_new_dir/my_artifact.csv", "nosuffix"])
@pytest.mark.parametrize("description", [None, "my description"])
def test_create_from_local_filepath(
    get_test_filepaths, key_is_virtual, key, description
):
    ln.settings.creation._artifact_use_virtual_keys = key_is_virtual
    is_in_registered_storage = get_test_filepaths[0]
    root_dir = get_test_filepaths[1]
    test_filepath = get_test_filepaths[3]
    suffix = get_test_filepaths[4]
    # this tests if insufficient information is being provided
    if key is None and not is_in_registered_storage and description is None:
        # this can fail because ln.context.track() might set a global run context
        # in that case, the File would have a run that's not None and the
        # error below wouldn't be thrown
        with pytest.raises(ValueError) as error:
            artifact = ln.Artifact(test_filepath, key=key, description=description)
        assert (
            error.exconly()
            == "ValueError: Pass one of key, run or description as a parameter"
        )
        return None
    elif key is not None and is_in_registered_storage:
        inferred_key = get_relative_path_to_directory(
            path=test_filepath, directory=root_dir
        ).as_posix()
        with pytest.raises(ValueError) as error:
            artifact = ln.Artifact(test_filepath, key=key, description=description)
        assert (
            error.exconly()
            == f"ValueError: The path '{test_filepath}' is already in registered"
            " storage"
            f" '{root_dir.resolve().as_posix()}' with key '{inferred_key}'\nYou"
            f" passed conflicting key '{key}': please move the file before"
            " registering it."
        )
        return None
    else:
        artifact = ln.Artifact(test_filepath, key=key, description=description)
        assert artifact._state.adding  # make sure that this is a new file in the db
    assert (
        artifact.description is None
        if description is None
        else artifact.description == description
    )
    assert artifact.suffix == suffix
    assert artifact.n_objects is None
    artifact.save()
    assert artifact.path.exists()

    if key is None:
        assert (
            artifact.key == f"my_dir/my_file{suffix}"
            if is_in_registered_storage
            else artifact.key is None
        )
        if is_in_registered_storage:
            assert artifact.storage.root == root_dir.resolve().as_posix()
            assert artifact.path == test_filepath.resolve()
        else:
            assert artifact.storage.root == lamindb_setup.settings.storage.root_as_str
            assert (
                artifact.path
                == lamindb_setup.settings.storage.root
                / f".lamindb/{artifact.uid}{suffix}"
            )
    else:
        assert artifact.key == key
        assert artifact._key_is_virtual == key_is_virtual
        if is_in_registered_storage:
            # this would only hit if the key matches the correct key
            assert artifact.storage.root == root_dir.resolve().as_posix()
            assert (
                artifact.path == root_dir / f"{key}{suffix}" == test_filepath.resolve()
            )
        else:
            # file is moved into default storage
            if key_is_virtual:
                assert (
                    artifact.path
                    == lamindb_setup.settings.storage.root
                    / f".lamindb/{artifact.uid}{suffix}"
                )
            else:
                assert artifact.path == lamindb_setup.settings.storage.root / key

    # only delete from storage if a file copy took place
    delete_from_storage = str(test_filepath.resolve()) != str(artifact.path)
    artifact.delete(permanent=True, storage=delete_from_storage)
    ln.settings.creation._artifact_use_virtual_keys = True


ERROR_MESSAGE = """\
ValueError: Currently don't support tracking folders outside one of the storage roots:
"""


@pytest.mark.parametrize("key", [None, "my_new_folder"])
def test_from_dir_many_artifacts(get_test_filepaths, key):
    is_in_registered_storage = get_test_filepaths[0]
    test_dirpath = get_test_filepaths[2]
    # the directory contains 3 files, two of them are duplicated
    if key is not None and is_in_registered_storage:
        with pytest.raises(ValueError) as error:
            ln.Artifact.from_dir(test_dirpath, key=key)
        assert error.exconly().startswith(
            "ValueError: The path"  # The path {data} is already in registered storage
        )
        return None
    else:
        artifacts = ln.Artifact.from_dir(test_dirpath, key=key)
    # we only return the duplicated ones
    hashes = [artifact.hash for artifact in artifacts if artifact.hash is not None]
    uids = [artifact.uid for artifact in artifacts]
    assert len(set(hashes)) == len(hashes)
    ln.UPath(test_dirpath).view_tree()
    # now save
    ln.save(artifacts)
    # now run again, because now we'll have hash-based lookup!
    artifacts = ln.Artifact.from_dir(test_dirpath, key=key)
    assert len(artifacts) == 2
    assert len(set(artifacts)) == len(hashes)
    queried_artifacts = ln.Artifact.filter(uid__in=uids).all()
    for artifact in queried_artifacts:
        artifact.delete(permanent=True, storage=False)


def test_delete_artifact(df):
    artifact = ln.Artifact.from_df(df, description="My test file to delete")
    artifact.save()
    assert artifact.visibility == 1
    assert artifact.key is None or artifact._key_is_virtual
    storage_path = artifact.path
    # trash behavior
    artifact.delete()
    assert storage_path.exists()
    assert ln.Artifact.filter(description="My test file to delete").first() is None
    assert ln.Artifact.filter(
        description="My test file to delete", visibility=-1
    ).first()
    # permanent delete
    artifact.delete(permanent=True)
    assert (
        ln.Artifact.filter(
            description="My test file to delete", visibility=None
        ).first()
        is None
    )
    assert not storage_path.exists()  # deletes from storage is key_is_virtual

    # test deleting artifact from non-managed storage
    artifact = ln.Artifact(
        "s3://lamindb-dev-datasets/file-to-test-for-delete.csv",
        description="My test file to delete from non-default storage",
    )
    artifact.save()
    assert artifact.storage.instance_uid != ln.setup.settings.instance.uid
    assert artifact.key is not None
    filepath = artifact.path
    with pytest.raises(IntegrityError) as e:
        artifact.delete()
    assert e.exconly().startswith(
        "lamindb.core.exceptions.IntegrityError: Cannot simply delete artifacts"
    )
    artifact.delete(storage=False, permanent=True)
    assert (
        ln.Artifact.filter(
            description="My test file to delete from non-default storage",
            visibility=None,
        ).first()
        is None
    )
    assert filepath.exists()


def test_delete_storage():
    with pytest.raises(FileNotFoundError):
        delete_storage(ln.settings.storage.root / "test-delete-storage")


# why does this run so long? in particular the first time?
@pytest.mark.parametrize(
    "filepath_str",
    ["s3://lamindb-ci/test-data/test.parquet", "s3://lamindb-ci/test-data/test.csv"],
)
@pytest.mark.parametrize("skip_check_exists", [False, True])
@pytest.mark.parametrize("skip_size_and_hash", [False, True])
def test_create_small_file_from_remote_path(
    filepath_str, skip_check_exists, skip_size_and_hash
):
    ln.settings.creation.artifact_skip_size_hash = skip_size_and_hash
    artifact = ln.Artifact(
        filepath_str,
        skip_check_exists=skip_check_exists,
    )
    artifact.save()
    # test cache()
    file_from_local = ln.Artifact(artifact.cache(), description="test")
    # test hash equivalency when computed on local machine
    if not skip_size_and_hash:
        assert file_from_local.hash == artifact.hash
        assert file_from_local._hash_type == "md5"
        assert artifact._hash_type == "md5"
    assert artifact.path.as_posix() == filepath_str
    assert artifact.load().iloc[0].tolist() == [
        0,
        "Abingdon island giant tortoise",
        "Chelonoidis abingdonii",
        106734,
        "ASM359739v1",
        "GCA_003597395.1",
        "Full genebuild",
        "-",
        "-",
    ]
    artifact.delete(permanent=True, storage=False)
    ln.settings.creation.artifact_skip_size_hash = False


def test_create_big_file_from_remote_path():
    previous_storage = ln.setup.settings.storage.root_as_str
    ln.settings.storage = "s3://lamindb-test"
    filepath_str = "s3://lamindb-test/human_immune.h5ad"
    artifact = ln.Artifact(filepath_str)
    assert artifact.key == "human_immune.h5ad"
    assert artifact._hash_type == "md5-2"
    ln.settings.storage = previous_storage


# -------------------------------------------------------------------------------------
# Storage
# -------------------------------------------------------------------------------------


@pytest.mark.parametrize("suffix", [".txt", "", None])
def test_auto_storage_key_from_artifact_uid(suffix):
    test_id = "abo389f"
    if suffix is None:
        with pytest.raises(AssertionError):
            auto_storage_key_from_artifact_uid(test_id, suffix, False)
    else:
        assert AUTO_KEY_PREFIX == ".lamindb/"
        storage_key = auto_storage_key_from_artifact_uid(test_id, suffix, False)
        assert storage_key == f"{AUTO_KEY_PREFIX}{test_id}{suffix}"


def test_storage_root_upath_equivalence():
    storage_root = UPath("s3://lamindb-ci")
    filepath = UPath("s3://lamindb-ci/test-data/Species.csv")
    assert filepath.parents[-1] == storage_root


def test_get_relative_path_to_directory():
    # upath on S3
    upath_root = UPath("s3://lamindb-ci")
    upath_directory1 = UPath("s3://lamindb-ci/test-data")  # no trailing slash
    upath_directory2 = UPath("s3://lamindb-ci/test-data/")  # trailing slash
    upath_file = UPath("s3://lamindb-ci/test-data/test.csv")
    assert (
        "test-data/test.csv"
        == get_relative_path_to_directory(upath_file, upath_root).as_posix()
    )
    assert (
        "test.csv"
        == get_relative_path_to_directory(upath_file, upath_directory1).as_posix()
    )
    assert (
        "test.csv"
        == get_relative_path_to_directory(upath_file, upath_directory2).as_posix()
    )
    # local path
    root = Path("/lamindb-ci")
    upath = Path("/lamindb-ci/test-data/test.csv")
    assert (
        "test-data/test.csv"
        == get_relative_path_to_directory(upath, directory=root).as_posix()
    )
    with pytest.raises(TypeError) as error:
        get_relative_path_to_directory(upath, directory=".")
    assert error.exconly() == "TypeError: Directory not of type Path or UPath"


def test_check_path_is_child_of_root():
    # UPath
    root = UPath("s3://lamindb-ci")
    upath = UPath("s3://lamindb-ci/test-data/test.csv")
    assert check_path_is_child_of_root(upath, root=root)
    upath2 = UPath("s3://lamindb-setup/test-data/test.csv")
    assert not check_path_is_child_of_root(upath2, root=root)
    # local path
    root = Path("/lamindb-ci")
    path = Path("/lamindb-ci/test-data/test.csv")
    assert check_path_is_child_of_root(path, root=root)
    path = Path("/lamindb-other/test-data/test.csv")
    assert not check_path_is_child_of_root(path, root=root)
    # Local & UPath
    root = UPath("s3://lamindb-ci")
    path = Path("/lamindb-ci/test-data/test.csv")
    assert not check_path_is_child_of_root(path, root=root)
    # different storage_options
    upath = UPath("s3://lamindb-ci/test-data/test.csv", cache_regions=True)
    assert upath.storage_options != root.storage_options
    assert check_path_is_child_of_root(upath, root=root)
    # the second level
    root = UPath("s3://lamindb-ci/test-data/")
    upath = UPath("s3://lamindb-ci/test-data/test/test.csv")
    assert check_path_is_child_of_root(upath, root=root)
    upath2 = UPath("s3://lamindb-ci/test-data-1/test/test.csv")
    assert not check_path_is_child_of_root(upath2, root=root)


def test_serialize_paths():
    fp_str = ln.core.datasets.anndata_file_pbmc68k_test().as_posix()
    fp_path = Path(fp_str)

    up_str = "s3://lamindb-ci/test-data/test.csv"
    up_upath = UPath(up_str)

    default_storage = settings._storage_settings.record
    using_key = None

    _, filepath, _, _, _ = process_data(
        "id", fp_str, None, None, default_storage, using_key, skip_existence_check=True
    )
    assert isinstance(filepath, LocalPathClasses)
    _, filepath, _, _, _ = process_data(
        "id", fp_path, None, None, default_storage, using_key, skip_existence_check=True
    )
    assert isinstance(filepath, LocalPathClasses)

    _, filepath, _, _, _ = process_data(
        "id", up_str, None, None, default_storage, using_key, skip_existence_check=True
    )
    assert isinstance(filepath, CloudPath)
    _, filepath, _, _, _ = process_data(
        "id",
        up_upath,
        None,
        None,
        default_storage,
        using_key,
        skip_existence_check=True,
    )
    assert isinstance(filepath, CloudPath)


def test_load_to_memory(tsv_file, zip_file, fcs_file):
    # tsv
    df = read_tsv(tsv_file)
    assert isinstance(df, pd.DataFrame)
    # fcs
    adata = read_fcs(fcs_file)
    assert isinstance(adata, ad.AnnData)
    # none
    load_to_memory(zip_file)

    with pytest.raises(TypeError) as error:
        ln.Artifact(True)
    assert error.exconly() == "TypeError: data has to be a string, Path, UPath"


def test_describe():
    ln.core.datasets.file_mini_csv()
    artifact = ln.Artifact("mini.csv", description="test")
    artifact.describe()


def test_zarr_upload_cache(adata):
    previous_storage = ln.setup.settings.storage.root_as_str
    ln.settings.storage = "s3://lamindb-test"

    def callback(*args, **kwargs):
        pass

    zarr_path = Path("./test_adata.zarr")
    write_adata_zarr(adata, zarr_path, callback)

    artifact = ln.Artifact(zarr_path, key="test_adata.zarr")
    assert artifact._accessor == "AnnData"
    assert artifact.n_objects >= 1
    artifact.save()

    assert isinstance(artifact.path, CloudPath)
    assert artifact.path.exists()
    assert zarr_is_adata(artifact.path)

    shutil.rmtree(artifact.cache())

    cache_path = settings._storage_settings.cloud_to_local_no_update(artifact.path)
    assert isinstance(artifact.load(), ad.AnnData)
    assert cache_path.is_dir()

    shutil.rmtree(cache_path)
    assert not cache_path.exists()
    artifact.cache()
    assert cache_path.is_dir()

    artifact.delete(permanent=True, storage=True)
    shutil.rmtree(zarr_path)

    # test zarr from memory
    artifact = ln.Artifact(adata, key="test_adata.anndata.zarr")
    assert artifact._local_filepath.is_dir()
    assert artifact._accessor == "AnnData"
    assert artifact.suffix == ".anndata.zarr"
    assert artifact.n_objects >= 1

    artifact.save()
    assert isinstance(artifact.path, CloudPath)
    assert artifact.path.exists()
    cache_path = settings._storage_settings.cloud_to_local_no_update(artifact.path)
    assert cache_path.is_dir()

    shutil.rmtree(cache_path)
    assert not cache_path.exists()

    artifact._memory_rep = None

    assert isinstance(artifact.load(), ad.AnnData)
    assert cache_path.is_dir()

    artifact.delete(permanent=True, storage=True)

    ln.settings.storage = previous_storage


def test_df_suffix(df):
    artifact = ln.Artifact.from_df(df, key="test_.parquet")
    assert artifact.suffix == ".parquet"

    with pytest.raises(ValueError) as error:
        artifact = ln.Artifact.from_df(df, key="test_.def")
    assert (
        error.exconly().partition(",")[0]
        == "ValueError: The suffix '.def' of the provided key is incorrect"
    )


def test_adata_suffix(adata):
    artifact = ln.Artifact.from_anndata(adata, key="test_.h5ad")
    assert artifact.suffix == ".h5ad"
    artifact = ln.Artifact.from_anndata(adata, format="h5ad", key="test_.h5ad")
    assert artifact.suffix == ".h5ad"
    artifact = ln.Artifact.from_anndata(adata, key="test_.zarr")
    assert artifact.suffix == ".zarr"

    with pytest.raises(ValueError) as error:
        artifact = ln.Artifact.from_anndata(adata, key="test_.def")
    assert (
        error.exconly().partition(",")[0]
        == "ValueError: Error when specifying AnnData storage format"
    )

    with pytest.raises(ValueError) as error:
        artifact = ln.Artifact.from_anndata(adata, key="test_")
    assert (
        error.exconly().partition(",")[0]
        == "ValueError: The suffix '' of the provided key is incorrect"
    )


def test_bulk_delete():
    report_path = Path("report.html")
    with open(report_path, "w") as f:
        f.write("a")
    _source_code_artifact_path = Path("code.py")
    with open(_source_code_artifact_path, "w") as f:
        f.write("b")
    environment_path = Path("environment.txt")
    with open(environment_path, "w") as f:
        f.write("c")
    report = ln.Artifact(report_path, description="Report")
    report.save()
    report_path.unlink()
    report_path = report.path
    _source_code_artifact = ln.Artifact(
        _source_code_artifact_path, description="Source"
    )
    _source_code_artifact.save()
    _source_code_artifact_path.unlink()
    _source_code_artifact_path = _source_code_artifact.path
    environment = ln.Artifact(environment_path, description="requirement.txt")
    environment.save()
    environment_path.unlink()
    environment_path = environment.path

    ln.Artifact.filter(
        id__in=[_source_code_artifact.id, environment.id, report.id]
    ).delete()

    assert (
        len(
            ln.Artifact.filter(
                id__in=[_source_code_artifact.id, environment.id, report.id]
            ).all()
        )
        == 0
    )
    assert (
        len(
            ln.Artifact.filter(
                id__in=[_source_code_artifact.id, environment.id, report.id],
                visibility=-1,
            )
            .distinct()
            .all()
        )
        == 3
    )

    ln.Artifact.filter(
        id__in=[_source_code_artifact.id, environment.id, report.id], visibility=-1
    ).delete(permanent=True)
    assert (
        len(
            ln.Artifact.filter(
                id__in=[_source_code_artifact.id, environment.id, report.id],
                visibility=None,
            )
            .distinct()
            .all()
        )
        == 0
    )

    assert not report_path.exists()
    assert not _source_code_artifact_path.exists()
    assert not environment_path.exists()
