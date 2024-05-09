import lamindb as ln
import pytest
from lamindb._save import prepare_error_message, store_artifacts


def test_bulk_save_and_update():
    label_names = [f"ULabel {i} new" for i in range(3)]
    labels = [ln.ULabel(name=name) for name in label_names]
    # test bulk creation of new records
    ln.save(labels)
    assert len(ln.ULabel.filter(name__in=label_names).distinct().all()) == 3
    labels[0].name = "ULabel 0 updated"
    # test bulk update of existing records
    ln.save(labels)
    assert len(ln.ULabel.filter(name__in=label_names).distinct().all()) == 2
    assert ln.ULabel.filter(name="ULabel 0 updated").one()


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


def test_save_data_object():
    import anndata as ad

    ln.core.datasets.file_mini_csv()
    artifact = ln.Artifact("mini.csv", description="test")
    artifact.save()
    assert artifact.path.exists()
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
