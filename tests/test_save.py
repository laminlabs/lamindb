import pytest

import lamindb as ln
from lamindb._save import prepare_error_message, store_files


def test_prepare_error_message():
    ln.dev.datasets.file_mini_csv()
    file = ln.File("mini.csv", description="test")
    exception = Exception("exception")

    error = prepare_error_message([], [file], exception)
    assert error.startswith(
        "The following entries have been successfully uploaded and committed to the database"  # noqa
    )

    error = prepare_error_message([file], [], exception)
    assert error.startswith("No entries were uploaded or committed to the database")


def test_zarr_upload_data_object():
    import anndata as ad

    ln.dev.datasets.file_mini_csv()
    file = ln.File("mini.csv", description="test")

    file.suffix = ".zarr"
    file._memory_rep = ad.AnnData()
    file.save()
    file.delete(permanent=True, storage=True)


def test_store_files_acid():
    ln.dev.datasets.file_mini_csv()
    file = ln.File("mini.csv", description="test")
    file._clear_storagekey = "test.csv"

    with pytest.raises(RuntimeError) as error:
        store_files([file])

    assert str(error.exconly()).startswith(
        "RuntimeError: No entries were uploaded or committed to the database."
    )
