import lamindb as ln
import pytest
from lamindb._save import prepare_error_message, store_artifacts


def test_prepare_error_message():
    ln.core.datasets.file_mini_csv()
    artifact = ln.Artifact("mini.csv", description="test")
    exception = Exception("exception")

    error = prepare_error_message([], [artifact], exception)
    assert error.startswith(
        "The following entries have been successfully uploaded and committed to the database"
    )

    error = prepare_error_message([artifact], [], exception)
    assert error.startswith("No entries were uploaded or committed to the database")


def test_zarr_upload_data_object():
    import anndata as ad

    ln.core.datasets.file_mini_csv()
    artifact = ln.Artifact("mini.csv", description="test")

    artifact.suffix = ".zarr"
    artifact._memory_rep = ad.AnnData()
    artifact.save()
    artifact.delete(permanent=True, storage=True)


def test_store_artifacts_acid():
    ln.core.datasets.file_mini_csv()
    artifact = ln.Artifact("mini.csv", description="test")
    artifact._clear_storagekey = "test.csv"

    with pytest.raises(RuntimeError) as error:
        store_artifacts([artifact], using_key=None)

    assert str(error.exconly()).startswith(
        "RuntimeError: No entries were uploaded or committed to the database."
    )
