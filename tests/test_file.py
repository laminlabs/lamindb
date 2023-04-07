from upath import UPath


def test_storage_root_upath_equivalence():
    storage_root = UPath("s3://lamindb-ci")
    filepath = UPath("s3://lamindb-ci/test-data/Species.csv")
    assert filepath.parents[-1] == storage_root
