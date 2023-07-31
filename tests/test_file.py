import shutil
from inspect import signature
from pathlib import Path

import anndata as ad
import lnschema_bionty as lb
import numpy as np
import pandas as pd
import pytest
from lamindb_setup.dev.upath import UPath

import lamindb as ln
from lamindb import _file
from lamindb._file import check_path_in_default_storage, get_relative_path_to_root

# how do we properly abstract out the default storage variable?
# currently, we're only mocking it through `default_storage` as
# set in conftest.py

lb.settings.species = "human"

df = pd.DataFrame({"feat1": [1, 2], "feat2": [3, 4]})

adata = ad.AnnData(
    X=np.array([[1, 2, 3], [4, 5, 6]]),
    obs=dict(Obs=["A", "B"]),
    var=dict(Feat=["MYC1", "TCF7", "GATA1"]),
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
    file.delete(storage=True)


def test_create_from_anndata_in_memory():
    file = ln.File.from_anndata(adata, var_ref=lb.Gene.symbol)
    assert file.accessor == "AnnData"
    assert hasattr(file, "_local_filepath")
    file.save()
    # check that the local filepath has been cleared
    assert not hasattr(file, "_local_filepath")
    feature_sets_queried = file.feature_sets.all()
    feature_list_queried = ln.Feature.filter(
        feature_sets__in=feature_sets_queried
    ).list()
    feature_list_queried = [feature.name for feature in feature_list_queried]
    assert set(feature_list_queried) == set(adata.obs.columns)
    feature_list_queried = lb.Gene.filter(feature_sets__in=feature_sets_queried).list()
    feature_list_queried = [feature.symbol for feature in feature_list_queried]
    assert set(feature_list_queried) == set(adata.var.index)
    file.delete(storage=True)


@pytest.mark.parametrize(
    "data", [adata, "s3://lamindb-test/scrnaseq_pbmc68k_tiny.h5ad"]
)
def test_create_from_anndata_in_storage(data):
    if isinstance(data, ad.AnnData):
        filepath = Path("test_adata.h5ad")
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
    scope="module", params=[(True, "./default_storage/"), (False, "./outside_storage/")]
)
def get_test_filepaths(request):  # -> Tuple[bool, Path, Path]
    isin_default_storage: bool = request.param[0]
    root_dir: Path = Path(request.param[1])
    test_dir = root_dir / "my_dir/"
    test_dir.mkdir(parents=True)
    test_filepath = test_dir / "my_file.csv"
    test_filepath.write_text("a")
    # return a boolean indicating whether test filepath is in default storage
    # and the test filepath
    yield (isin_default_storage, test_dir, test_filepath)
    shutil.rmtree(test_dir)


# this tests the basic (non-provenance-related) metadata
@pytest.mark.parametrize("key", [None, "my_new_dir/my_file.csv"])
@pytest.mark.parametrize("name", [None, "my name"])
def test_create_from_local_filepath(get_test_filepaths, key, name):
    isin_default_storage = get_test_filepaths[0]
    test_filepath = get_test_filepaths[2]
    file = ln.File(test_filepath, key=key, name=name)
    assert file.description is None if name is None else file.description == name
    assert file.suffix == ".csv"
    if key is None:
        assert (
            file.key == "my_dir/my_file.csv"
            if isin_default_storage
            else file.key is None
        )
    else:
        assert file.key == key
    assert file.storage.root == Path("./default_storage").resolve().as_posix()
    assert file.hash == "DMF1ucDxtqgxw5niaXcmYQ"
    if isin_default_storage and key is None:
        assert str(test_filepath.resolve()) == str(file.path())

    # need to figure out how to save!
    # now save it
    # file.save()
    # assert file.path().exists()
    # file.delete(storage=True)


def test_local_path_load():
    local_filepath = Path("tests/test-files/pbmc68k.h5ad").resolve()

    file = ln.File(local_filepath)
    assert local_filepath == file._local_filepath
    assert local_filepath == file.path()
    assert local_filepath == file.stage()

    adata = ad.read(local_filepath)
    file = ln.File(adata)
    assert file._memory_rep is adata
    assert file.load() is adata
    assert file._local_filepath.resolve() == file.stage() == file.path()


def test_init_from_directory(get_test_filepaths):
    isin_default_storage = get_test_filepaths[0]
    test_dirpath = get_test_filepaths[1]
    if isin_default_storage:
        storage_root = None
    else:
        storage_root = test_dirpath.parent
    records = ln.File.from_dir(test_dirpath, storage_root=storage_root)
    assert len(records) == 1
    # also execute tree
    if isin_default_storage:
        ln.File.tree(test_dirpath.name)


def test_delete(get_test_filepaths):
    test_filepath = get_test_filepaths[2]
    file = ln.File(test_filepath, name="My test file to delete")
    file.save()
    storage_path = file.path()
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
    assert file.path().as_posix() == filepath_str
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


def test_get_name_suffix_from_filepath():
    # based on https://stackoverflow.com/questions/31890341/clean-way-to-get-the-true-stem-of-a-path-object  # noqa
    dataset = [
        ("a", "a", ""),
        ("a.txt", "a", ".txt"),
        ("archive.tar.gz", "archive", ".tar.gz"),
        ("directory/file", "file", ""),
        ("d.x.y.z/f.a.b.c", "f", ".a.b.c"),
        ("logs/date.log.txt", "date", ".log.txt"),
    ]
    for path, _, suffix in dataset:
        filepath = Path(path)
        assert suffix == "".join(filepath.suffixes)


def test_storage_root_upath_equivalence():
    storage_root = UPath("s3://lamindb-ci")
    filepath = UPath("s3://lamindb-ci/test-data/Species.csv")
    assert filepath.parents[-1] == storage_root


def test_get_relative_path_to_root():
    # upath on S3
    root = UPath("s3://lamindb-ci")
    upath = UPath("s3://lamindb-ci/test-data/test.csv")
    assert (
        "test-data/test.csv" == get_relative_path_to_root(upath, root=root).as_posix()
    )
    # local path
    root = Path("/lamindb-ci")
    upath = Path("/lamindb-ci/test-data/test.csv")
    assert (
        "test-data/test.csv" == get_relative_path_to_root(upath, root=root).as_posix()
    )


def test_check_path_in_default_storage():
    # UPath
    root = UPath("s3://lamindb-ci")
    upath = UPath("s3://lamindb-ci/test-data/test.csv")
    assert check_path_in_default_storage(upath, root=root)
    upath2 = UPath("s3://lamindb-setup/test-data/test.csv")
    assert not check_path_in_default_storage(upath2, root=root)
    # local path
    root = Path("/lamindb-ci")
    path = Path("/lamindb-ci/test-data/test.csv")
    assert check_path_in_default_storage(path, root=root)
    path = Path("/lamindb-other/test-data/test.csv")
    assert not check_path_in_default_storage(path, root=root)
    # Local & UPath
    root = UPath("s3://lamindb-ci")
    path = Path("/lamindb-ci/test-data/test.csv")
    assert not check_path_in_default_storage(path, root=root)
