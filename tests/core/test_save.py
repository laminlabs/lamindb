import lamindb as ln
import pytest
from lamindb.models.save import prepare_error_message, store_artifacts


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
    assert ln.ULabel.get(name="ULabel 0 updated")


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
    ln.core.datasets.file_mini_csv()
    artifact = ln.Artifact("mini.csv", description="test")
    artifact.save()
    assert artifact.path.exists()
    artifact.delete(permanent=True, storage=True)


def test_store_artifacts_acid():
    ln.core.datasets.file_mini_csv()
    artifact = ln.Artifact("mini.csv", description="test")
    artifact._clear_storagekey = "test.csv"
    # errors on check_and_attempt_clearing
    with pytest.raises(RuntimeError):
        artifact.save()

    with pytest.raises(RuntimeError) as error:
        store_artifacts([artifact], using_key=None)
    assert str(error.exconly()).startswith(
        "RuntimeError: The following entries have been successfully uploaded"
    )

    artifact.delete(permanent=True)


def test_save_parents():
    import bionty as bt

    records = bt.CellLine.from_values(["HEPG2", "HUVEC"])
    ln.save(records)
    assert bt.CellLine.get("4ea731nb").parents.df().shape[0] == 1
    bt.CellLine.filter().delete()
