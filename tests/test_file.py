import shutil
from inspect import signature
from pathlib import Path

import anndata as ad
import lamindb_setup
import lnschema_bionty as lb
import numpy as np
import pandas as pd
import pytest
from django.db.models.deletion import ProtectedError
from lamindb_setup.dev.upath import (
    CloudPath,
    LocalPathClasses,
    UPath,
    extract_suffix_from_path,
)

import lamindb as ln
from lamindb import _file
from lamindb._file import (
    check_path_is_child_of_root,
    get_hash,
    get_relative_path_to_directory,
    process_data,
)
from lamindb.dev.storage._zarr import write_adata_zarr
from lamindb.dev.storage.file import (
    AUTO_KEY_PREFIX,
    auto_storage_key_from_id_suffix,
    delete_storage,
    load_to_memory,
    read_fcs,
    read_tsv,
)

# how do we properly abstract out the default storage variable?
# currently, we're only mocking it through `default_storage` as
# set in conftest.py

ln.settings.verbosity = "success"
lb.settings.organism = "human"

df = pd.DataFrame({"feat1": [1, 2], "feat2": [3, 4]})

adata = ad.AnnData(
    X=np.array([[1, 2, 3], [4, 5, 6]]),
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
    class_methods = ["view_tree", "from_dir", "from_df", "from_anndata"]
    for name in class_methods:
        setattr(Mock, name, getattr(_file, name))
        assert signature(getattr(Mock, name)) == _file.SIGS.pop(name)
    # methods
    for name, sig in _file.SIGS.items():
        assert signature(getattr(_file, name)) == sig


@pytest.fixture(
    scope="module",
    params=[
        # tuple of isin_existing_storage, path, suffix
        (True, "./default_storage/", ".csv"),
        (True, "./default_storage/", ""),
        (True, "./registered_storage/", ".csv"),
        (True, "./registered_storage/", ""),
        (False, "./nonregistered_storage/", ".csv"),
        (False, "./nonregistered_storage/", ""),
    ],
)
def get_test_filepaths(request):  # -> Tuple[bool, Path, Path, Path, str]
    import lamindb as ln

    isin_existing_storage: bool = request.param[0]
    root_dir: Path = Path(request.param[1])
    suffix: str = request.param[2]
    if isin_existing_storage:
        # ensure that it's actually registered
        if ln.Storage.filter(root=root_dir.resolve().as_posix()).one_or_none() is None:
            ln.Storage(root=root_dir.resolve().as_posix(), type="local").save()
    else:
        assert (
            ln.Storage.filter(root=root_dir.resolve().as_posix()).one_or_none() is None
        )
    test_dir = root_dir / "my_dir/"
    test_dir.mkdir(parents=True)
    test_filepath = test_dir / f"my_file{suffix}"
    test_filepath.write_text(str(test_filepath))
    # create a duplicated file
    test_filepath1 = test_dir / f"my_file1{suffix}"
    test_filepath1.write_text(str(test_filepath))
    # create a non-duplicated file
    test_filepath2 = test_dir / f"my_file2{suffix}"
    test_filepath2.write_text(str(test_filepath2))
    # return a boolean indicating whether test filepath is in default storage
    # and the test filepath
    yield (isin_existing_storage, root_dir, test_dir, test_filepath, suffix)
    shutil.rmtree(test_dir)


def test_is_new_version_of_versioned_file():
    # attempt to create a file with an invalid version
    with pytest.raises(ValueError) as error:
        file = ln.File(df, description="test", version=0)
    assert (
        error.exconly()
        == "ValueError: `version` parameter must be `None` or `str`, e.g., '0.1', '1',"
        " '2', etc."
    )

    # create a versioned file
    file = ln.File(df, description="test", version="1")
    assert file.version == "1"

    assert file.path.exists()  # because of cache file already exists
    file.save()
    assert file.path.exists()

    with pytest.raises(ValueError) as error:
        file_v2 = ln.File(adata, is_new_version_of=file, version="1")
    assert error.exconly() == "ValueError: Please increment the previous version: '1'"

    # create new file from old file
    file_v2 = ln.File(adata, is_new_version_of=file)
    assert file.version == "1"
    assert file.initial_version_id is None  # initial file has initial_version_id None
    assert file_v2.initial_version_id == file.id
    assert file_v2.version == "2"
    assert file_v2.key is None
    assert file_v2.description == "test"

    file_v2.save()
    assert file_v2.path.exists()

    # create new file from newly versioned file
    df.iloc[0, 0] = 0
    file_v3 = ln.File(df, description="test1", is_new_version_of=file_v2)
    assert file_v3.initial_version_id == file.id
    assert file_v3.version == "3"
    assert file_v3.description == "test1"

    with pytest.raises(TypeError) as error:
        ln.File(df, description="test1a", is_new_version_of=ln.Transform())
    assert error.exconly() == "TypeError: is_new_version_of has to be of type ln.File"

    # test that reference file cannot be deleted
    with pytest.raises(ProtectedError):
        file.delete(permanent=True, storage=True)
    file_v2.delete(permanent=True, storage=True)
    file.delete(permanent=True, storage=True)

    # extra kwargs
    with pytest.raises(ValueError):
        ln.File(df, description="test1b", extra_kwarg="extra")

    # > 1 args
    with pytest.raises(ValueError) as error:
        ln.File(df, df)
    assert error.exconly() == "ValueError: Only one non-keyword arg allowed: data"

    # AUTO_KEY_PREFIX
    with pytest.raises(ValueError) as error:
        ln.File(df, key=".lamindb/test_df.parquet")
    assert error.exconly() == "ValueError: Key cannot start with .lamindb/"


def test_is_new_version_of_unversioned_file():
    # unversioned file
    file = ln.File(df, description="test2")
    assert file.initial_version_id is None
    assert file.version is None

    # what happens if we don't save the old file?
    # add a test for it!
    file.save()

    # create new file from old file
    new_file = ln.File(adata, is_new_version_of=file)
    assert file.version == "1"
    assert file.initial_version is None
    assert new_file.initial_version_id == file.id
    assert new_file.version == "2"
    assert new_file.description == file.description

    file.delete(permanent=True, storage=True)


# also test legacy name parameter (got removed by description)
def test_create_from_dataframe():
    file = ln.File(df, description="test1")
    assert file.description == "test1"
    assert file.key is None
    assert file.accessor == "DataFrame"
    assert hasattr(file, "_local_filepath")
    file.save()
    # can't do backed
    with pytest.raises(ValueError):
        file.backed()
    # check that the local filepath has been cleared
    assert not hasattr(file, "_local_filepath")
    file.delete(permanent=True, storage=True)


def test_create_from_dataframe_using_from_df():
    description = "my description"
    file = ln.File.from_df(df, key="folder/hello.parquet", description=description)
    assert file._feature_sets == {}
    with pytest.raises(ValueError):
        file.features["columns"]
    ln.save(ln.Feature.from_df(df))
    file = ln.File.from_df(df, description=description)
    file = ln.File.from_df(df, key="folder/hello.parquet", description=description)
    # mere access test right now
    file.features["columns"]
    assert file.description == description
    assert file.accessor == "DataFrame"
    assert hasattr(file, "_local_filepath")
    assert file.key == "folder/hello.parquet"
    assert file.key_is_virtual
    assert file.uid in file.path.as_posix()
    file.save()
    # check that the local filepath has been cleared
    assert not hasattr(file, "_local_filepath")
    feature_set_queried = file.feature_sets.get()  # exactly one
    feature_list_queried = ln.Feature.filter(feature_sets=feature_set_queried).list()
    feature_list_queried = [feature.name for feature in feature_list_queried]
    assert set(feature_list_queried) == set(df.columns)
    feature_set_queried.delete()
    file.delete(permanent=True, storage=True)
    ln.Feature.filter(name__in=["feat1", "feat2"]).delete()


def test_create_from_anndata_in_memory():
    ln.save(lb.Gene.from_values(adata.var.index, "symbol"))
    ln.save(ln.Feature.from_df(adata.obs))
    file = ln.File.from_anndata(adata, description="test", field=lb.Gene.symbol)
    assert file.accessor == "AnnData"
    assert hasattr(file, "_local_filepath")
    file.save()
    # check that the local filepath has been cleared
    assert not hasattr(file, "_local_filepath")
    feature_sets_queried = file.feature_sets.all()
    features_queried = ln.Feature.filter(feature_sets__in=feature_sets_queried).all()
    assert set(features_queried.list("name")) == set(adata.obs.columns)
    genes_queried = lb.Gene.filter(feature_sets__in=feature_sets_queried).all()
    assert set(genes_queried.list("symbol")) == set(adata.var.index)
    feature_sets_queried.delete()
    features_queried.delete()
    genes_queried.delete()
    file.delete(permanent=True, storage=True)


@pytest.mark.parametrize(
    "data", [adata, "s3://lamindb-test/scrnaseq_pbmc68k_tiny.h5ad"]
)
def test_create_from_anndata_in_storage(data):
    if isinstance(data, ad.AnnData):
        filepath = Path("./default_storage/test_adata.h5ad")
        data.write(filepath)
    else:
        previous_storage = ln.setup.settings.storage.root_as_str
        ln.settings.storage = "s3://lamindb-test"
        filepath = data
    file = ln.File.from_anndata(filepath, field=lb.Gene.symbol)
    assert file.accessor == "AnnData"
    assert hasattr(file, "_local_filepath")
    file.save()
    # check that the local filepath has been cleared
    assert not hasattr(file, "_local_filepath")
    if isinstance(data, ad.AnnData):
        filepath.unlink()
    else:
        ln.settings.storage = previous_storage


# this tests the basic (non-provenance-related) metadata
@pytest.mark.parametrize("key_is_virtual", [True, False])
@pytest.mark.parametrize("key", [None, "my_new_dir/my_file.csv", "nosuffix"])
@pytest.mark.parametrize("description", [None, "my description"])
def test_create_from_local_filepath(
    get_test_filepaths, key_is_virtual, key, description
):
    ln.settings.file_use_virtual_keys = key_is_virtual
    isin_existing_storage = get_test_filepaths[0]
    root_dir = get_test_filepaths[1]
    test_filepath = get_test_filepaths[3]
    suffix = get_test_filepaths[4]
    # this tests if insufficient information is being provided
    if key is None and not isin_existing_storage and description is None:
        # this can fail because ln.track() might set a global run context
        # in that case, the File would have a run that's not None and the
        # error below wouldn't be thrown
        with pytest.raises(ValueError) as error:
            file = ln.File(test_filepath, key=key, description=description)
        assert (
            error.exconly()
            == "ValueError: Pass one of key, run or description as a parameter"
        )
        return None
    elif key is not None and isin_existing_storage:
        inferred_key = get_relative_path_to_directory(
            path=test_filepath, directory=root_dir
        ).as_posix()
        with pytest.raises(ValueError) as error:
            file = ln.File(test_filepath, key=key, description=description)
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
        file = ln.File(test_filepath, key=key, description=description)
        assert file._state.adding  # make sure that this is a new file in the db
    assert (
        file.description is None
        if description is None
        else file.description == description
    )
    assert file.suffix == suffix

    file.save()
    assert file.path.exists()

    if key is None:
        assert (
            file.key == f"my_dir/my_file{suffix}"
            if isin_existing_storage
            else file.key is None
        )
        if isin_existing_storage:
            assert file.storage.root == root_dir.resolve().as_posix()
            assert file.path == test_filepath.resolve()
        else:
            assert file.storage.root == lamindb_setup.settings.storage.root_as_str
            assert (
                file.path
                == lamindb_setup.settings.storage.root / f".lamindb/{file.uid}{suffix}"
            )
    else:
        assert file.key == key
        assert file.key_is_virtual == key_is_virtual
        if isin_existing_storage:
            # this would only hit if the key matches the correct key
            assert file.storage.root == root_dir.resolve().as_posix()
            assert file.path == root_dir / f"{key}{suffix}" == test_filepath.resolve()
        else:
            # file is moved into default storage
            if key_is_virtual:
                assert (
                    file.path
                    == lamindb_setup.settings.storage.root
                    / f".lamindb/{file.uid}{suffix}"
                )
            else:
                assert file.path == lamindb_setup.settings.storage.root / key

    # only delete from storage if a file copy took place
    delete_from_storage = str(test_filepath.resolve()) != str(file.path)
    file.delete(permanent=True, storage=delete_from_storage)
    ln.settings.file_use_virtual_keys = True


def test_local_path_load():
    local_filepath = ln.dev.datasets.anndata_file_pbmc68k_test().resolve()

    file = ln.File(local_filepath, description="test")
    assert local_filepath == file._local_filepath
    assert local_filepath == file.path
    assert local_filepath == file.stage()

    adata = ad.read(local_filepath)
    file = ln.File(adata, description="test")
    assert file._memory_rep is adata
    assert file.load() is adata
    assert file._local_filepath.resolve() == file.stage() == file.path


ERROR_MESSAGE = """\
ValueError: Currently don't support tracking folders outside one of the storage roots:
"""


@pytest.mark.parametrize("key", [None, "my_new_folder"])
def test_from_dir(get_test_filepaths, key):
    isin_existing_storage = get_test_filepaths[0]
    # root_dir = get_test_filepaths[1]
    test_dirpath = get_test_filepaths[2]
    # the directory contains 3 files, two of them are duplicated
    if key is not None and isin_existing_storage:
        with pytest.raises(ValueError) as error:
            ln.File.from_dir(test_dirpath, key=key)
        assert error.exconly().startswith(
            "ValueError: The path"  # The path {data} is already in registered storage
        )
        return None
    else:
        files = ln.File.from_dir(test_dirpath, key=key)
    # we only return the duplicated ones
    hashes = [file.hash for file in files if file.hash is not None]
    assert len(set(hashes)) == len(hashes)
    ln.UPath(test_dirpath).view_tree()
    # now save
    ln.save(files)
    ln.File.view_tree()
    ln.File.filter().all().view_tree()
    # now run again, because now we'll have hash-based lookup!
    files = ln.File.from_dir(test_dirpath, key=key)
    assert len(files) == 2
    assert len(set(hashes)) == len(hashes)
    for file in files:
        file.delete(permanent=True, storage=False)


def test_delete(get_test_filepaths):
    test_filepath = get_test_filepaths[3]
    file = ln.File(test_filepath, description="My test file to delete")
    file.save()
    storage_path = file.path
    file.delete(permanent=True, storage=True)
    assert ln.File.filter(description="My test file to delete").first() is None
    assert not Path(storage_path).exists()


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
    ln.settings.upon_file_create_skip_size_hash = skip_size_and_hash
    file = ln.File(
        filepath_str,
        skip_check_exists=skip_check_exists,
    )
    file.save()
    # test stage()
    file_from_local = ln.File(file.stage(), description="test")
    # test hash equivalency when computed on local machine
    if not skip_size_and_hash:
        assert file_from_local.hash == file.hash
        assert file_from_local.hash_type == "md5"
        assert file.hash_type == "md5"
    assert file.path.as_posix() == filepath_str
    assert file.load().iloc[0].tolist() == [
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
    file.delete(permanent=True, storage=False)
    ln.settings.upon_file_create_skip_size_hash = False


def test_create_big_file_from_remote_path():
    previous_storage = ln.setup.settings.storage.root_as_str
    ln.settings.storage = "s3://lamindb-test"
    filepath_str = "s3://lamindb-test/human_immune.h5ad"
    file = ln.File(filepath_str)
    assert file.key == "human_immune.h5ad"
    assert file.hash.endswith("-2")
    assert file.hash_type == "md5-n"
    ln.settings.storage = previous_storage


# -------------------------------------------------------------------------------------
# Storage
# -------------------------------------------------------------------------------------


def test_extract_suffix_from_path():
    # this is a dataset of path, stem, suffix tuples
    dataset = [
        ("a", "a", ""),
        ("a.txt", "a", ".txt"),
        ("archive.tar.gz", "archive", ".tar.gz"),
        ("directory/file", "file", ""),
        ("d.x.y.z/f.b.c", "f", ".c"),
        ("d.x.y.z/f.a.b.c", "f", ".c"),
        ("logs/date.log.txt", "date", ".txt"),
        ("salmon.merged.gene_counts.tsv", "salmon.merged.gene_counts", ".tsv"),
        ("salmon.merged.gene_counts.tsv.gz", "salmon.merged.gene_counts", ".tsv.gz"),
    ]
    for path, _, suffix in dataset:
        filepath = Path(path)
        assert suffix == extract_suffix_from_path(filepath)


@pytest.mark.parametrize("suffix", [".txt", "", None])
def test_auto_storage_key_from_id_suffix(suffix):
    test_id = "abo389f"
    if suffix is None:
        with pytest.raises(AssertionError):
            auto_storage_key_from_id_suffix(test_id, suffix)
    else:
        assert AUTO_KEY_PREFIX == ".lamindb/"
        storage_key = auto_storage_key_from_id_suffix(test_id, suffix)
        storage_key == f"{AUTO_KEY_PREFIX}{test_id}{suffix}"


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
    # default root
    assert not check_path_is_child_of_root(path)


def test_serialize_paths():
    fp_str = ln.dev.datasets.anndata_file_pbmc68k_test().as_posix()
    fp_path = Path(fp_str)

    up_str = "s3://lamindb-ci/test-data/test.csv"
    up_upath = UPath(up_str)

    _, filepath, _, _, _ = process_data(
        "id", fp_str, None, None, skip_existence_check=True
    )
    assert isinstance(filepath, LocalPathClasses)
    _, filepath, _, _, _ = process_data(
        "id", fp_path, None, None, skip_existence_check=True
    )
    assert isinstance(filepath, LocalPathClasses)

    _, filepath, _, _, _ = process_data(
        "id", up_str, None, None, skip_existence_check=True
    )
    assert isinstance(filepath, CloudPath)
    _, filepath, _, _, _ = process_data(
        "id", up_upath, None, None, skip_existence_check=True
    )
    assert isinstance(filepath, CloudPath)


def test_load_to_memory():
    # tsv
    pd.DataFrame([1, 2]).to_csv("test.tsv", sep="\t")
    df = read_tsv("test.tsv")
    assert isinstance(df, pd.DataFrame)
    # fcs
    adata = read_fcs(ln.dev.datasets.file_fcs())
    assert isinstance(adata, ad.AnnData)
    # other
    pd.DataFrame([1, 2]).to_csv("test.zrad", sep="\t")
    with pytest.raises(NotADirectoryError):
        load_to_memory("test.zrad")
    # none
    pd.DataFrame([1, 2]).to_csv("test.zip", sep="\t")
    load_to_memory("test.zip")
    assert get_hash("test.zip", suffix=".zip", check_hash=False)[0] is None
    UPath("test.tsv").unlink()
    UPath("test.zrad").unlink()
    UPath("test.zip").unlink()

    with pytest.raises(NotImplementedError) as error:
        ln.File(True)
    assert (
        error.exconly()
        == "NotImplementedError: Do not know how to create a file object from True,"
        " pass a filepath instead!"
    )


def test_delete_storage():
    with pytest.raises(FileNotFoundError):
        delete_storage(UPath("test"))


def test_describe():
    ln.dev.datasets.file_mini_csv()
    file = ln.File("mini.csv", description="test")
    file.describe()


def test_file_zarr():
    with open("test.zarr", "w") as f:
        f.write("zarr")
    assert get_hash(UPath("test.zarr"), suffix=".zarr") is None
    file = ln.File("test.zarr", description="test-zarr")
    with pytest.raises(RuntimeError) as error:
        file.stage()
    assert (
        error.exconly()
        == "RuntimeError: zarr object can't be staged, please use load() or stream()"
    )  # noqa
    file.save()
    file.delete(permanent=True, storage=False)
    UPath("test.zarr").unlink()


def test_zarr_folder_upload():
    previous_storage = ln.setup.settings.storage.root_as_str
    ln.settings.storage = "s3://lamindb-test"

    def callback(*args, **kwargs):
        pass

    zarr_path = Path("./test_adata.zrad")
    write_adata_zarr(adata, zarr_path, callback)

    file = ln.File(zarr_path, key="test_adata.zrad")
    file.save()

    assert isinstance(file.path, CloudPath) and file.path.exists()

    file.delete(permanent=True, storage=True)
    delete_storage(zarr_path)
    ln.settings.storage = previous_storage


def test_df_suffix():
    file = ln.File(df, key="test_.parquet")
    assert file.suffix == ".parquet"

    with pytest.raises(ValueError) as error:
        file = ln.File(df, key="test_.def")
    assert (
        error.exconly().partition(",")[0]
        == "ValueError: The suffix '.def' of the provided key is incorrect"
    )


def test_adata_suffix():
    file = ln.File(adata, key="test_.h5ad")
    assert file.suffix == ".h5ad"
    file = ln.File(adata, format="h5ad", key="test_.h5ad")
    assert file.suffix == ".h5ad"
    file = ln.File(adata, key="test_.zarr")
    assert file.suffix == ".zarr"
    file = ln.File(adata, key="test_.zrad")
    assert file.suffix == ".zrad"
    file = ln.File(adata, format="zrad", key="test_.zrad")
    assert file.suffix == ".zrad"

    with pytest.raises(ValueError) as error:
        file = ln.File(adata, key="test_.def")
    assert (
        error.exconly().partition(",")[0]
        == "ValueError: Error when specifying AnnData storage format"
    )

    with pytest.raises(ValueError) as error:
        file = ln.File(adata, format="h5ad", key="test.zrad")
    assert (
        error.exconly().partition(",")[0]
        == "ValueError: The suffix '.zrad' of the provided key is incorrect"
    )

    with pytest.raises(ValueError) as error:
        file = ln.File(adata, key="test_")
    assert (
        error.exconly().partition(",")[0]
        == "ValueError: The suffix '' of the provided key is incorrect"
    )
