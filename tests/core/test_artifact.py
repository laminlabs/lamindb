"""Artifact tests.

Also see `test_artifact_folders.py` for tests of folder-like artifacts.

"""

# ruff: noqa: F811

import shutil
import sys
from pathlib import Path, PurePosixPath
from types import ModuleType

import anndata as ad
import bionty as bt
import lamindb as ln
import lamindb_setup
import mudata as md
import numpy as np
import pandas as pd
import pytest
import yaml  # type: ignore
from _dataset_fixtures import (  # noqa
    get_small_adata,
    get_small_mdata,
    get_small_sdata,
    get_small_soma_experiment,
)
from lamindb.core._settings import settings
from lamindb.core.loaders import load_fcs, load_to_memory, load_tsv
from lamindb.core.storage.paths import (
    AUTO_KEY_PREFIX,
    auto_storage_key_from_artifact_uid,
    delete_storage,
)
from lamindb.errors import (
    FieldValidationError,
    InvalidArgument,
)
from lamindb.models.artifact import (
    check_path_is_child_of_root,
    data_is_scversedatastructure,
    data_is_soma_experiment,
    get_relative_path_to_directory,
    process_data,
)
from lamindb_setup.core.upath import (
    CloudPath,
    LocalPathClasses,
    UPath,
    extract_suffix_from_path,
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
    if request.param == "get_small_adata":
        return request.getfixturevalue("get_small_adata")
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
def yaml_file():
    filepath = Path("test.yaml")
    dct = {"a": 1, "b": 2}
    with open(filepath, "w") as f:
        yaml.dump(dct, f)
    yield filepath
    filepath.unlink()


@pytest.fixture(scope="module")
def fcs_file():
    return ln.core.datasets.file_fcs_alpert19()


@pytest.fixture(scope="module")
def mudata_file(get_small_mdata):
    filepath = Path("test.h5mu")
    get_small_mdata.write(filepath)
    yield filepath
    filepath.unlink()


@pytest.fixture(scope="module")
def spatialdata_file(get_small_sdata):
    filepath = Path("test.zarr")
    get_small_sdata.write(filepath)
    yield filepath
    shutil.rmtree(filepath)


@pytest.fixture(scope="function")
def soma_experiment_file(get_small_soma_experiment, clean_soma_files):
    yield "test.tiledbsoma"


def test_data_is_anndata_paths():
    assert data_is_scversedatastructure("something.h5ad", "AnnData")
    assert data_is_scversedatastructure("something.anndata.zarr", "AnnData")
    assert data_is_scversedatastructure(
        "s3://somewhere/something.anndata.zarr", "AnnData"
    )
    assert not data_is_scversedatastructure("s3://somewhere/something.zarr", "AnnData")


def test_data_is_mudata_paths():
    assert data_is_scversedatastructure("something.h5mu", "MuData")
    assert data_is_scversedatastructure("something.mudata.zarr", "MuData")


def test_data_is_spatialdata_paths():
    assert data_is_scversedatastructure("something.spatialdata.zarr", "SpatialData")


def test_data_is_soma_experiment_paths():
    assert data_is_soma_experiment("something.tiledbsoma")


@pytest.mark.parametrize(
    "data,data_type,expected",
    [
        ("get_small_adata", "AnnData", True),
        ("get_small_mdata", "MuData", True),
        ("get_small_sdata", "SpatialData", True),
        ("get_small_adata", "MuData", False),
        ("get_small_mdata", "AnnData", False),
        ("get_small_sdata", "AnnData", False),
        ("get_small_adata", None, True),
        (pd.DataFrame(), "AnnData", False),
        (None, "AnnData", False),
        (None, None, False),
    ],
)
def test_data_is_scversedatastructure(request, data, data_type, expected):
    if isinstance(data, str) and data.startswith("get_small_"):
        data = request.getfixturevalue(data)

    assert data_is_scversedatastructure(data, data_type) == expected


def test_data_is_soma_experiment(get_small_soma_experiment, clean_soma_files):
    assert data_is_soma_experiment(get_small_soma_experiment)


def test_basic_validation():
    # extra kwargs
    with pytest.raises(FieldValidationError):
        ln.Artifact("testpath.csv", description="test1b", extra_kwarg="extra")

    # > 1 args
    with pytest.raises(ValueError) as error:
        ln.Artifact("testpath.csv", "testpath.csv")
    assert error.exconly() == "ValueError: Only one non-keyword arg allowed: data"

    # AUTO_KEY_PREFIX in key
    with pytest.raises(ValueError) as error:
        ln.Artifact.from_df(df, key=".lamindb/test_df.parquet")
    assert (
        error.exconly()
        == f"ValueError: Do not pass key that contains a managed storage path in `{AUTO_KEY_PREFIX}`"
    )

    # path that contains AUTO_KEY_PREFIX
    with pytest.raises(ValueError) as error:
        ln.Artifact(".lamindb/test_df.parquet", description="Test")
    assert (
        error.exconly()
        == f"ValueError: Do not pass path inside the `{AUTO_KEY_PREFIX}` directory."
    )


def test_revise_artifact(df):
    # attempt to create a file with an invalid version
    with pytest.raises(ValueError) as error:
        artifact = ln.Artifact.from_df(df, description="test", version=0)
    assert (
        error.exconly()
        == "ValueError: `version` parameter must be `None` or `str`, e.g., '0.1', '1',"
        " '2', etc."
    )

    # create a file and tag it with a version
    key = "my-test-dataset.parquet"
    artifact = ln.Artifact.from_df(df, key=key, description="test", version="1")
    assert artifact.version == "1"
    assert artifact.uid.endswith("0000")
    assert artifact.path.exists()  # because of cache file already exists
    artifact.save()
    assert artifact.path.exists()
    assert artifact.suffix == ".parquet"

    with pytest.raises(ValueError) as error:
        artifact_r2 = ln.Artifact.from_df(df, revises=artifact, version="1")
    assert error.exconly() == "ValueError: Please increment the previous version: '1'"

    # create new file from old file
    df.iloc[0, 0] = 99  # mutate dataframe so that hash lookup doesn't trigger
    artifact_r2 = ln.Artifact.from_df(df, revises=artifact)
    assert artifact_r2.stem_uid == artifact.stem_uid
    assert artifact_r2.uid.endswith("0001")
    # call this again
    artifact_r2 = ln.Artifact.from_df(df, revises=artifact)
    assert artifact_r2.uid.endswith("0001")
    assert artifact_r2.stem_uid == artifact.stem_uid
    assert artifact_r2.version is None
    assert artifact_r2.key == key
    assert artifact.suffix == ".parquet"
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
    df.iloc[0, 0] = 100  # mutate dataframe so that hash lookup doesn't trigger
    artifact_r3 = ln.Artifact.from_df(df, description="test1", key=key, version="2")
    assert artifact_r3.uid.endswith("0002")
    assert artifact_r3.stem_uid == artifact.stem_uid
    assert artifact_r3.key == key
    assert artifact_r3.version == "2"
    assert artifact_r3.description == "test1"
    assert artifact_r3.is_latest
    assert artifact_r2.is_latest
    artifact_r3.save()
    # now r2 is no longer the latest version, but need to re-fresh from db
    artifact_r2 = ln.Artifact.get(artifact_r2.uid)
    assert not artifact_r2.is_latest

    # re-create based on hash while providing a different key
    artifact_new = ln.Artifact.from_df(
        df,
        description="test1 updated",
        key="my-test-dataset1.parquet",
    )
    assert artifact_new == artifact_r3
    assert artifact_new.key == key  # old key
    assert artifact_new.description == "test1 updated"

    with pytest.raises(TypeError) as error:
        ln.Artifact.from_df(df, description="test1a", revises=ln.ULabel(name="test"))
    assert error.exconly() == "TypeError: `revises` has to be of type `Artifact`"

    artifact_r3.delete(permanent=True)
    artifact_r2.delete(permanent=True)
    artifact.delete(permanent=True)

    # unversioned file
    artifact = ln.Artifact.from_df(df, description="test2")
    assert artifact.version is None

    # what happens if we don't save the old file?
    # add a test for it!
    artifact.save()

    # create new file from old file
    df.iloc[0, 0] = 101  # mutate dataframe so that hash lookup doesn't trigger
    new_artifact = ln.Artifact.from_df(df, revises=artifact)
    assert artifact.version is None
    assert new_artifact.stem_uid == artifact.stem_uid
    assert new_artifact.version is None
    assert new_artifact.description == artifact.description

    artifact.delete(permanent=True)


def test_create_from_dataframe(df):
    artifact = ln.Artifact.from_df(df, description="test1")
    assert artifact.description == "test1"
    assert artifact.key is None
    assert artifact.otype == "DataFrame"
    assert artifact.kind == "dataset"
    assert artifact.n_observations == 2
    assert hasattr(artifact, "_local_filepath")
    artifact.key = "my-test-dataset"  # try changing key
    with pytest.raises(InvalidArgument) as error:
        artifact.save()
    assert (
        error.exconly()
        == "lamindb.errors.InvalidArgument: The suffix '' of the provided key is incorrect, it should be '.parquet'."
    )
    artifact.key = None  # restore
    artifact.suffix = ".whatever"  # try changing suffix
    with pytest.raises(InvalidArgument) as error:
        artifact.save()
    assert (
        error.exconly()
        == "lamindb.errors.InvalidArgument: Changing the `.suffix` of an artifact is not allowed! You tried to change it from '.parquet' to '.whatever'."
    )
    artifact.suffix = ".parquet"
    artifact.save()
    # check that the local filepath has been cleared
    assert not hasattr(artifact, "_local_filepath")
    del artifact

    # now get an artifact from the database
    artifact = ln.Artifact.get(description="test1")

    artifact.suffix = ".whatever"  # try changing suffix
    with pytest.raises(InvalidArgument) as error:
        artifact.save()
    assert (
        error.exconly()
        == "lamindb.errors.InvalidArgument: Changing the `.suffix` of an artifact is not allowed! You tried to change it from '.parquet' to '.whatever'."
    )
    artifact.suffix = ".parquet"

    # coming from `key is None` that setting a key with different suffix is not allowed
    artifact.key = "my-test-dataset.suffix"
    with pytest.raises(InvalidArgument) as error:
        artifact.save()
    assert (
        error.exconly()
        == "lamindb.errors.InvalidArgument: The suffix '.suffix' of the provided key is incorrect, it should be '.parquet'."
    )

    # coming from `key is None` test with no suffix
    artifact.key = "my-test-dataset"
    with pytest.raises(InvalidArgument) as error:
        artifact.save()
    assert (
        error.exconly()
        == "lamindb.errors.InvalidArgument: The suffix '' of the provided key is incorrect, it should be '.parquet'."
    )

    # try a joint update where both suffix and key are changed, this previously enabled to create a corrupted state
    artifact.key = "my-test-dataset"
    artifact.suffix = ""
    with pytest.raises(InvalidArgument) as error:
        artifact.save()
    assert (
        error.exconly()
        == "lamindb.errors.InvalidArgument: Changing the `.suffix` of an artifact is not allowed! You tried to change it from '.parquet' to ''."
    )

    # because this is a parquet artifact, we can set a key with a .parquet suffix
    artifact.suffix = ".parquet"  # restore proper suffix
    artifact.key = "my-test-dataset.parquet"
    artifact.save()
    assert artifact.key == "my-test-dataset.parquet"

    # coming from a .parquet key, test changing the key to no suffix
    artifact.key = "my-test-dataset"
    with pytest.raises(InvalidArgument) as error:
        artifact.save()
    assert (
        error.exconly()
        == "lamindb.errors.InvalidArgument: The suffix '' of the provided key is incorrect, it should be '.parquet'."
    )

    artifact.delete(permanent=True)


def test_create_from_anndata(get_small_adata, adata_file, df):
    with pytest.raises(ValueError) as error:
        ln.Artifact.from_anndata(df, description="test1")
    assert (
        "data has to be an AnnData object or a path to AnnData-like" in error.exconly()
    )

    for i, _a in enumerate([get_small_adata, adata_file]):
        artifact = ln.Artifact.from_anndata(_a, description="test1")
        assert artifact.description == "test1"
        assert artifact.key is None
        assert artifact.otype == "AnnData"
        assert artifact.kind == "dataset"
        assert artifact.n_observations == 2
        if i == 0:
            assert hasattr(artifact, "_local_filepath")
            artifact.save()
            # check that the local filepath has been cleared
            assert not hasattr(artifact, "_local_filepath")
            artifact.delete(permanent=True)


def test_create_from_mudata(get_small_mdata, mudata_file, adata_file):
    with pytest.raises(ValueError) as error:
        ln.Artifact.from_mudata(adata_file, description="test1")
    assert "data has to be a MuData object or a path to MuData-like" in error.exconly()

    for m in [get_small_mdata, mudata_file]:
        af = ln.Artifact.from_mudata(m, description="test1")
        assert af.description == "test1"
        assert af.key is None
        assert af.otype == "MuData"
        assert af.kind == "dataset"
        if isinstance(m, md.MuData):
            assert af.n_observations == 2


def test_create_from_spatialdata(
    get_small_sdata, spatialdata_file, adata_file, ccaplog
):
    with pytest.raises(ValueError) as error:
        ln.Artifact.from_spatialdata(adata_file, description="test1")
    assert (
        "data has to be a SpatialData object or a path to SpatialData-like"
        in error.exconly()
    )

    for s in [get_small_sdata, spatialdata_file]:
        af = ln.Artifact(s, description="test1")
        assert af.description == "test1"
        assert af.key is None
        assert af.otype == "SpatialData"
        assert af.kind is None
        # n_observations not defined
    assert "data is a SpatialData, please use .from_spatialdata()" in ccaplog.text
    for s in [get_small_sdata, spatialdata_file]:
        af = ln.Artifact.from_spatialdata(s, description="test1")
        assert af.description == "test1"
        assert af.key is None
        assert af.otype == "SpatialData"
        assert af.kind == "dataset"
        # n_observations not defined


def test_create_from_soma_experiment(
    soma_experiment_file, clean_soma_files, adata_file
):
    with pytest.raises(ValueError) as error:
        ln.Artifact.from_tiledbsoma(adata_file, description="test1")
    assert (
        "data has to be a SOMA Experiment object or a path to SOMA Experiment store."
        in error.exconly()
    )

    af = ln.Artifact.from_tiledbsoma(soma_experiment_file, description="test1")
    assert af.description == "test1"
    assert af.key is None
    assert af.otype == "tiledbsoma"
    assert af.n_observations == 3


@pytest.mark.parametrize(
    "data",
    ["get_small_adata"],
    indirect=True,
)
def test_create_from_anndata_in_storage(data):
    artifact = ln.Artifact.from_anndata(
        data, description="test_create_from_anndata_memory"
    )
    assert artifact.n_observations == data.n_obs
    assert artifact.otype == "AnnData"
    assert hasattr(artifact, "_local_filepath")
    artifact.save()
    # check that the local filepath has been cleared
    assert not hasattr(artifact, "_local_filepath")


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
    suffix = get_test_filepaths[4]  # path suffix
    if key is not None:
        key_suffix = extract_suffix_from_path(
            PurePosixPath(key), arg_name="key"
        )  # key suffix
    else:
        key_suffix = None
    # this tests if insufficient information is being provided
    if key is None and not is_in_registered_storage and description is None:
        # this can fail because ln.track() might set a global run context
        # in that case, the Artifact would have a run that's not None and the
        # error below wouldn't be thrown
        with pytest.raises(ValueError) as error:
            artifact = ln.Artifact(test_filepath, key=key, description=description)
        assert (
            error.exconly()
            == "ValueError: Pass one of key, run or description as a parameter"
        )
        return None
    elif key is not None and suffix != key_suffix:
        with pytest.raises(InvalidArgument) as error:
            artifact = ln.Artifact(test_filepath, key=key, description=description)
        assert error.exconly() == (
            f"lamindb.errors.InvalidArgument: The passed path's suffix '{suffix}' must match the passed key's suffix '{key_suffix}'."
        )
        return None
    elif key is not None and is_in_registered_storage:
        inferred_key = get_relative_path_to_directory(
            path=test_filepath, directory=root_dir
        ).as_posix()
        try:
            artifact = ln.Artifact(test_filepath, key=key, description=description)
        except InvalidArgument as error:
            assert str(error) == (
                f"The path '{test_filepath}' is already in registered storage"
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
    assert artifact.n_files is None
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
        # changing non-virtual key is not allowed
        if not key_is_virtual:
            with pytest.raises(InvalidArgument):
                artifact.key = "new_key"
                artifact.save()
            # need to change the key back to the original key
            artifact.key = key
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


@pytest.mark.parametrize("key", [None, "my_new_folder"])
def test_from_dir_many_artifacts(get_test_filepaths, key):
    is_in_registered_storage = get_test_filepaths[0]
    test_dirpath = get_test_filepaths[2]
    # the directory contains 3 files, two of them are duplicated
    if key is not None and is_in_registered_storage:
        with pytest.raises(InvalidArgument) as error:
            ln.Artifact.from_dir(test_dirpath, key=key)
        assert error.exconly().startswith(
            "lamindb.errors.InvalidArgument: The path"  # The path {data} is already in registered storage
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


def test_delete_and_restore_artifact(df):
    artifact = ln.Artifact.from_df(df, description="My test file to delete").save()
    assert artifact.branch_id == 1
    assert artifact.key is None or artifact._key_is_virtual
    storage_path = artifact.path
    # trash behavior
    artifact.delete()
    assert storage_path.exists()
    assert artifact.branch_id == -1
    assert ln.Artifact.filter(description="My test file to delete").first() is None
    assert ln.Artifact.filter(
        description="My test file to delete", branch_id=-1
    ).first()
    # implicit restore from trash
    artifact_restored = ln.Artifact.from_df(df, description="My test file to delete")
    assert artifact_restored.branch_id == 1
    assert artifact_restored == artifact
    # permanent delete
    artifact.delete(permanent=True)
    assert (
        ln.Artifact.filter(description="My test file to delete", branch_id=None).first()
        is None
    )
    assert not storage_path.exists()  # deletes from storage is key_is_virtual


def test_delete_storage():
    with pytest.raises(FileNotFoundError):
        delete_storage(ln.settings.storage.root / "test-delete-storage")


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
    # str
    root = "s3://lamindb-ci"
    upath = "s3://lamindb-ci/test-data/test.csv"
    assert check_path_is_child_of_root(upath, root=root)
    # str different protocols
    root = "prot1://lamindb-ci"
    upath = "prot2://lamindb-ci/test-data/test.csv"
    assert not check_path_is_child_of_root(upath, root=root)
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
    # http
    assert check_path_is_child_of_root(
        "https://raw.githubusercontent.com/laminlabs/lamindb/refs/heads/main/README.md",
        root="https://raw.githubusercontent.com",
    )
    # s3 with endpoint
    assert not check_path_is_child_of_root(
        "s3://bucket/key?endpoint_url=http://localhost:8000",
        root="s3://bucket/",
    )
    assert not check_path_is_child_of_root(
        "s3://bucket/key/",
        root="s3://bucket/?endpoint_url=http://localhost:8000",
    )
    assert check_path_is_child_of_root(
        "s3://bucket/key?endpoint_url=http://localhost:8000",
        root="s3://bucket?endpoint_url=http://localhost:8000",
    )
    assert check_path_is_child_of_root(
        UPath("s3://bucket/key", endpoint_url="http://localhost:8000"),
        root="s3://bucket?endpoint_url=http://localhost:8000",
    )


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


def test_load_to_memory(tsv_file, zip_file, fcs_file, yaml_file):
    # tsv
    df = load_tsv(tsv_file)
    assert isinstance(df, pd.DataFrame)
    # fcs
    adata = load_fcs(str(fcs_file))
    assert isinstance(adata, ad.AnnData)
    # error
    with pytest.raises(NotImplementedError):
        load_to_memory(zip_file)
    # check that it is a path
    assert isinstance(load_to_memory("./somefile.rds"), UPath)
    # yaml
    dct = load_to_memory(yaml_file)
    assert dct["a"] == 1
    assert dct["b"] == 2

    with pytest.raises(TypeError) as error:
        ln.Artifact(True)
    assert error.exconly() == "TypeError: data has to be a string, Path, UPath"


def test_describe():
    ln.core.datasets.file_mini_csv()
    artifact = ln.Artifact("mini.csv", description="test")
    artifact.describe()

    # test describing from a remote instance with less modules
    artifact = ln.Artifact.using("laminlabs/lamin-site-assets").filter().first()
    artifact.describe()


def test_df_suffix(df):
    artifact = ln.Artifact.from_df(df, key="test_.parquet")
    assert artifact.suffix == ".parquet"

    with pytest.raises(InvalidArgument) as error:
        artifact = ln.Artifact.from_df(df, key="test_.def")
    assert (
        error.exconly().partition(",")[0]
        == "lamindb.errors.InvalidArgument: The passed key's suffix '.def' must match the passed path's suffix '.parquet'."
    )


def test_adata_suffix(get_small_adata):
    artifact = ln.Artifact.from_anndata(get_small_adata, key="test_.h5ad")
    assert artifact.suffix == ".h5ad"
    artifact = ln.Artifact.from_anndata(
        get_small_adata, format="h5ad", key="test_.h5ad"
    )
    assert artifact.suffix == ".h5ad"
    artifact = ln.Artifact.from_anndata(get_small_adata, key="test_.zarr")
    assert artifact.suffix == ".zarr"

    with pytest.raises(ValueError) as error:
        artifact = ln.Artifact.from_anndata(get_small_adata, key="test_.def")
    assert (
        error.exconly().partition(",")[0]
        == "ValueError: Error when specifying AnnData storage format"
    )

    with pytest.raises(InvalidArgument) as error:
        artifact = ln.Artifact.from_anndata(get_small_adata, key="test_")
    assert (
        error.exconly().partition(",")[0]
        == "lamindb.errors.InvalidArgument: The passed key's suffix '' must match the passed path's suffix '.h5ad'."
    )


def test_bulk_delete():
    report_path = Path("report.html")
    report_path.write_text("a")
    environment_path = Path("environment.txt")
    environment_path.write_text("c")
    report = ln.Artifact(report_path, description="Report").save()
    report_path.unlink()
    report_path = report.path
    environment = ln.Artifact(environment_path, description="requirement.txt").save()
    environment_path.unlink()
    environment_path = environment.path

    ln.Artifact.filter(id__in=[environment.id, report.id]).delete()

    assert len(ln.Artifact.filter(id__in=[environment.id, report.id]).all()) == 0

    # the 2 artifacts are in trash now
    assert (
        len(
            ln.Artifact.filter(
                id__in=[environment.id, report.id],
                branch_id=-1,
            )
            .distinct()
            .all()
        )
        == 2
    )

    ln.Artifact.filter(id__in=[environment.id, report.id], branch_id=-1).delete(
        permanent=True
    )
    # now they're gone
    assert (
        len(
            ln.Artifact.filter(
                id__in=[environment.id, report.id],
                branch_id=None,
            )
            .distinct()
            .all()
        )
        == 0
    )

    assert not report_path.exists()
    assert not environment_path.exists()


@pytest.mark.parametrize("module_name", ["mudata", "spatialdata"])
def test_no_unnecessary_imports(df, module_name: str) -> None:
    if module_name in sys.modules:
        del sys.modules[module_name]

    af = ln.Artifact.from_df(df, description="to delete").save()

    loaded_packages = []
    for name, module in sys.modules.items():
        if isinstance(module, ModuleType) and not name.startswith("_"):
            if "." not in name:
                loaded_packages.append(name)

    assert module_name not in sorted(loaded_packages)

    # Cleanup and restore imports to ensure that other tests still run smoothly
    af.delete(permanent=True)
    import mudata  # noqa
    import spatialdata  # noqa


def test_artifact_get_tracking(df):
    artifact = ln.Artifact.from_df(df, key="df.parquet").save()

    transform = ln.Transform(key="test track artifact via get").save()
    run = ln.Run(transform).save()

    assert (
        ln.Artifact.get(key="df.parquet", is_run_input=run) in run.input_artifacts.all()
    )

    artifact.delete(permanent=True)
    transform.delete()
