# ruff: noqa: F811

import lamindb as ln
import pytest
from _dataset_fixtures import (  # noqa
    get_mini_csv,
)
from lamindb.models.save import prepare_error_message, store_artifacts


def test_bulk_save_and_update():
    label_names = [f"Record {i} new" for i in range(3)]
    labels = [ln.Record(name=name) for name in label_names]
    # test bulk creation of new records
    ln.save(labels)
    assert len(ln.Record.filter(name__in=label_names).distinct()) == 3
    labels[0].name = "Record 0 updated"
    # test bulk update of existing records
    ln.save(labels)
    assert len(ln.Record.filter(name__in=label_names).distinct()) == 2
    assert ln.Record.get(name="Record 0 updated")


def test_prepare_error_message(get_mini_csv):
    artifact = ln.Artifact(get_mini_csv, description="test")
    exception = Exception("exception")

    error = prepare_error_message([], [artifact], exception)
    assert error.startswith(
        "The following entries have been successfully uploaded and committed to the database"
    )

    error = prepare_error_message([artifact], [], exception)
    assert error.startswith("No entries were uploaded or committed to the database")


def test_save_data_object(get_mini_csv):
    artifact = ln.Artifact(get_mini_csv, description="test")
    artifact.save()
    assert artifact.path.exists()
    artifact.delete(permanent=True, storage=True)


def test_store_artifacts_acid(get_mini_csv):
    artifact = ln.Artifact(get_mini_csv, description="test")
    artifact._clear_storagekey = "test.csv"
    # errors on check_and_attempt_clearing
    with pytest.raises(FileNotFoundError):
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
    assert bt.CellLine.get("4ea731nb").parents.to_dataframe().shape[0] == 1
    bt.CellLine.filter().delete(permanent=True)


def test_save_batch_size():
    label_names = [f"Record {i} batch_size" for i in range(3)]
    labels = [ln.Record(name=name) for name in label_names]
    # test bulk creation of new records with batch size
    ln.save(labels, batch_size=2)
    assert ln.Record.filter(name__in=label_names).distinct().count() == 3


def test_bulk_resave_trashed_records():
    import bionty as bt

    # first create records from public source
    records = bt.Ethnicity.from_values(["asian", "white"]).save()
    # one parent is also created
    ethnicities = bt.Ethnicity.filter()
    assert ethnicities.count() == 3
    # soft delete the records including parent
    ethnicities.delete()
    # then create them again from public source
    # the new records will now have the same uids as they are hashed from the ontology_ids
    assert bt.Ethnicity.filter().count() == 0
    new_records = bt.Ethnicity.from_values(["asian", "white", "african"])
    assert new_records[0].branch_id == 1
    assert new_records[0].uid == records[0].uid
    # after saving, the trashed records should be restored
    new_records.save()
    assert new_records[0].branch_id == 1
    ethnicities = bt.Ethnicity.filter()
    # the parent should also be restored
    assert ethnicities.count() == 4

    # clean up
    ethnicities.delete(permanent=True)
