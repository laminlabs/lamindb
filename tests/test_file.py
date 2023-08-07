import shutil
from inspect import signature
from pathlib import Path

import anndata as ad
import lamindb_setup
import lnschema_bionty as lb
import numpy as np
import pandas as pd
import pytest
from lamindb_setup.dev.upath import UPath

import lamindb as ln
from lamindb import _file
from lamindb._file import (
    check_path_is_child_of_root,
    get_relative_path_to_directory,
    process_data,
)
from lamindb.dev.storage.file import extract_suffix_from_path

# how do we properly abstract out the default storage variable?
# currently, we're only mocking it through `default_storage` as
# set in conftest.py

ln.settings.verbosity = 3
lb.settings.species = "human"

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
    class_methods = ["tree", "from_dir", "from_df", "from_anndata"]
    for name in class_methods:
        setattr(Mock, name, getattr(_file, name))
        assert signature(getattr(Mock, name)) == _file.SIGS.pop(name)
    # methods
    for name, sig in _file.SIGS.items():
        assert signature(getattr(_file, name)) == sig


@pytest.mark.parametrize("name", [None, "my name"])
def test_create_from_dataframe(name):
    file = ln.File(df, name=name)
    assert file.description is None if name is None else file.description == name
    assert file.key is None
    assert file.accessor == "DataFrame"
    assert hasattr(file, "_local_filepath")
    file.save()
    # check that the local filepath has been cleared
    assert not hasattr(file, "_local_filepath")
    file.delete(storage=True)


@pytest.mark.parametrize("description", [None, "my name"])
def test_create_from_dataframe_using_from_df(description):
    file = ln.File.from_df(df, description=description)
    assert file._feature_sets == {}
    with pytest.raises(ValueError):
        file.features["columns"]
    ln.save(ln.Feature.from_df(df))
    file = ln.File.from_df(df, description=description)
    # mere access test right now
    file.features["columns"]
    assert file.description == description
    assert file.key is None
    assert file.accessor == "DataFrame"
    assert hasattr(file, "_local_filepath")
    file.save()
    # check that the local filepath has been cleared
    assert not hasattr(file, "_local_filepath")
    feature_set_queried = file.feature_sets.get()  # exactly one
    feature_list_queried = ln.Feature.filter(feature_sets=feature_set_queried).list()
    feature_list_queried = [feature.name for feature in feature_list_queried]
    assert set(feature_list_queried) == set(df.columns)
    feature_set_queried.delete()
    file.delete(storage=True)
    ln.Feature.filter(name__in=["feat1", "feat2"]).delete()


def test_create_from_anndata_in_memory():
    ln.save(ln.Feature.from_df(adata.obs))
    file = ln.File.from_anndata(adata, var_ref=lb.Gene.symbol)
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
    file.delete(storage=True)


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
    file = ln.File.from_anndata(filepath, var_ref=lb.Gene.symbol)
    assert file.accessor == "AnnData"
    assert hasattr(file, "_local_filepath")
    file.save()
    # check that the local filepath has been cleared
    assert not hasattr(file, "_local_filepath")
    if isinstance(data, ad.AnnData):
        filepath.unlink()
    else:
        ln.settings.storage = previous_storage


@pytest.fixture(
    scope="module",
    params=[
        (True, "./default_storage/"),
        (True, "./registered_storage/"),
        (False, "./nonregistered_storage/"),
    ],
)
def get_test_filepaths(request):  # -> Tuple[bool, Path, Path, Path]
    isin_existing_storage: bool = request.param[0]
    root_dir: Path = Path(request.param[1])
    if isin_existing_storage:
        # ensure that it's actually registered
        if ln.Storage.filter(root=root_dir.resolve().as_posix()).one_or_none() is None:
            ln.Storage(root=root_dir.resolve().as_posix(), type="local").save()
    test_dir = root_dir / "my_dir/"
    test_dir.mkdir(parents=True)
    test_filepath = test_dir / "my_file.csv"
    test_filepath.write_text("a")
    # return a boolean indicating whether test filepath is in default storage
    # and the test filepath
    yield (isin_existing_storage, root_dir, test_dir, test_filepath)
    shutil.rmtree(test_dir)


# this tests the basic (non-provenance-related) metadata
@pytest.mark.parametrize("key", [None, "my_new_dir/my_file.csv"])
@pytest.mark.parametrize("description", [None, "my description"])
def test_create_from_local_filepath(get_test_filepaths, key, description):
    isin_existing_storage = get_test_filepaths[0]
    root_dir = get_test_filepaths[1]
    test_filepath = get_test_filepaths[3]
    file = ln.File(test_filepath, key=key, description=description)
    assert (
        file.description is None
        if description is None
        else file.description == description
    )
    assert file.suffix == ".csv"
    if key is None:
        assert (
            file.key == "my_dir/my_file.csv"
            if isin_existing_storage
            else file.key is None
        )
        if isin_existing_storage:
            assert file.storage.root == root_dir.resolve().as_posix()
        else:
            assert file.storage.root == lamindb_setup.settings.storage.root_as_str
    else:
        assert file.key == key
        assert file.storage.root == lamindb_setup.settings.storage.root_as_str
    assert file.hash == "DMF1ucDxtqgxw5niaXcmYQ"

    # test that the file didn't move
    if isin_existing_storage and key is None:
        assert str(test_filepath.resolve()) == str(file.path)

    # now, save the file
    file.save()
    print(file.path)
    assert file.path.exists()

    # only delete from storage if a file copy took place
    delete_from_storage = str(test_filepath.resolve()) != str(file.path)
    file.delete(storage=delete_from_storage)


def test_local_path_load():
    local_filepath = ln.dev.datasets.anndata_file_pbmc68k_test().resolve()

    file = ln.File(local_filepath)
    assert local_filepath == file._local_filepath
    assert local_filepath == file.path
    assert local_filepath == file.stage()

    adata = ad.read(local_filepath)
    file = ln.File(adata)
    assert file._memory_rep is adata
    assert file.load() is adata
    assert file._local_filepath.resolve() == file.stage() == file.path


ERROR_MESSAGE = """\
ValueError: Currently don't support tracking folders outside one of the storage roots:
"""


@pytest.mark.parametrize("key", [None, "my_new_folder"])
def test_init_from_directory(get_test_filepaths, key):
    test_dirpath = get_test_filepaths[2]
    records = ln.File.from_dir(test_dirpath, key=key)
    assert len(records) == 1
    ln.File.tree(test_dirpath)


def test_delete(get_test_filepaths):
    test_filepath = get_test_filepaths[3]
    file = ln.File(test_filepath, description="My test file to delete")
    file.save()
    storage_path = file.path
    file.delete(storage=True)
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
    file_from_local = ln.File(file.stage())
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
    file.delete(storage=False)
    ln.settings.upon_file_create_skip_size_hash = False


def test_create_big_file_from_remote_path():
    previous_storage = ln.setup.settings.storage.root_as_str
    ln.settings.storage = "s3://lamindb-test"
    filepath_str = "s3://lamindb-test/human_immune.h5ad"
    file = ln.File(filepath_str)
    assert file.hash.endswith("-2")
    assert file.hash_type == "md5-n"
    ln.settings.storage = previous_storage


def test_inherit_relations():
    with open("test-inherit1", "w") as f:
        f.write("file1")
    with open("test-inherit2", "w") as f:
        f.write("file2")

    file1 = ln.File("test-inherit1")
    file1.save()
    file2 = ln.File("test-inherit2")
    file2.save()

    label_names = [f"Project {i}" for i in range(3)]
    labels = [ln.Label(name=label_name) for label_name in label_names]
    ln.save(labels)

    cell_line_names = [f"Cell line {i}" for i in range(3)]
    cell_lines = [lb.CellLine(name=name) for name in cell_line_names]
    ln.save(cell_lines)

    file2.labels.set(labels)
    file2.cell_lines.set(cell_lines)

    assert file1.labels.exists() is False
    file1.inherit_relations(file2, ["labels"])
    assert file1.labels.count() == file2.labels.count()
    assert file1.cell_lines.exists() is False
    file1.inherit_relations(file2)
    assert file1.cell_lines.count() == file2.cell_lines.count()

    with pytest.raises(KeyError):
        file1.inherit_relations(file2, ["not_exist_field"])

    for label in labels:
        label.delete()
    for cell_line in cell_lines:
        cell_line.delete()
    file1.delete(storage=True)
    file2.delete(storage=True)


# -------------------------------------------------------------------------------------
# Storage
# -------------------------------------------------------------------------------------


def test_extract_suffix_from_path():
    dataset = [
        ("a", "a", ""),
        ("a.txt", "a", ".txt"),
        ("archive.tar.gz", "archive", ".tar.gz"),
        ("directory/file", "file", ""),
        ("d.x.y.z/f.a.b.c", "f", ".c"),
        ("logs/date.log.txt", "date", ".log.txt"),
        ("salmon.merged.gene_counts.tsv", "salmon.merged.gene_counts", ".tsv"),
        ("salmon.merged.gene_counts.tsv.gz", "salmon.merged.gene_counts", ".tsv.gz"),
    ]
    for path, _, suffix in dataset:
        filepath = Path(path)
        assert suffix == extract_suffix_from_path(filepath)


def test_storage_root_upath_equivalence():
    storage_root = UPath("s3://lamindb-ci")
    filepath = UPath("s3://lamindb-ci/test-data/Species.csv")
    assert filepath.parents[-1] == storage_root


def test_get_relative_path_to_directory():
    # upath on S3
    root = UPath("s3://lamindb-ci")
    upath = UPath("s3://lamindb-ci/test-data/test.csv")
    assert (
        "test-data/test.csv"
        == get_relative_path_to_directory(upath, directory=root).as_posix()
    )
    # local path
    root = Path("/lamindb-ci")
    upath = Path("/lamindb-ci/test-data/test.csv")
    assert (
        "test-data/test.csv"
        == get_relative_path_to_directory(upath, directory=root).as_posix()
    )


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


def test_serialize_paths():
    fp_str = ln.dev.datasets.anndata_file_pbmc68k_test().as_posix()
    fp_path = Path(fp_str)

    up_str = "s3://lamindb-ci/test-data/test.csv"
    up_upath = UPath(up_str)

    _, filepath, _, _, _ = process_data("id", fp_str, None, skip_existence_check=True)
    assert isinstance(filepath, Path) and not isinstance(filepath, UPath)
    _, filepath, _, _, _ = process_data("id", fp_path, None, skip_existence_check=True)
    assert isinstance(filepath, Path) and not isinstance(filepath, UPath)

    _, filepath, _, _, _ = process_data("id", up_str, None, skip_existence_check=True)
    assert isinstance(filepath, UPath)
    _, filepath, _, _, _ = process_data("id", up_upath, None, skip_existence_check=True)
    assert isinstance(filepath, UPath)
