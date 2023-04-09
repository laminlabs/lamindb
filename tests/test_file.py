from pathlib import Path

from upath import UPath

from lamindb._file import get_check_path_in_storage, get_relative_path_to_root


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


def test_get_check_path_in_storage():
    # UPath
    root = UPath("s3://lamindb-ci")
    upath = UPath("s3://lamindb-ci/test-data/test.csv")
    assert get_check_path_in_storage(upath, root=root)
    upath2 = UPath("s3://lamindb-setup/test-data/test.csv")
    assert not get_check_path_in_storage(upath2, root=root)
    # local path
    root = Path("/lamindb-ci")
    path = Path("/lamindb-ci/test-data/test.csv")
    assert get_check_path_in_storage(path, root=root)
    path = Path("/lamindb-other/test-data/test.csv")
    assert not get_check_path_in_storage(path, root=root)
    # Local & UPath
    root = UPath("s3://lamindb-ci")
    path = Path("/lamindb-ci/test-data/test.csv")
    assert not get_check_path_in_storage(path, root=root)
