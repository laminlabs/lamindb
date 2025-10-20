import bionty as bt
import lamindb as ln
import pytest
from django.db import IntegrityError


@pytest.mark.parametrize("permanent", [True, False])
def test_delete_qs(permanent):
    """Test deletion behavior for small (1) and large (>=2) querysets.

    Small querysets delete individually, large ones trigger bulk delete."""
    ln.settings.creation.search_names = False
    labels = [ln.Record(name=f"label_{i}") for i in range(3)]
    ln.settings.creation.search_names = True
    ln.save(labels)
    ln.Record.filter(name__startswith="label_").delete(permanent=permanent)
    assert ln.Record.filter(name__startswith="label_", branch_id=-1).count() == (
        0 if permanent else 3
    )
    assert ln.ULabel.filter(name__startswith="label_").count() == 0


def test_delete_record():
    record_type = ln.Record(name="TestType", is_type=True).save()
    record = ln.Record(name="test_record", type=record_type).save()
    record2 = ln.Record(name="test_record", type=record_type).save()
    assert record == record2
    with pytest.raises(IntegrityError):
        # raise the unique constraint error when trying to create a duplicate record
        ln.Record(name="test_record", type=record_type, _skip_validation=True).save()
    record.delete()
    # because `record` is now in trash, we can create a new record with the same name and type
    record2 = ln.Record(name="test_record", type=record_type).save()
    assert record != record2
    record2.delete(permanent=True)
    record.delete(permanent=True)
    record_type.delete(permanent=True)


def test_recreate_soft_deleted_record():
    # testing soft delete and recreate with postgres (sqlite is tested in curators/test_records.py)
    # soft delete a record, then recreate it with some changes
    record = bt.Ethnicity.from_source(ontology_id="HANCESTRO:0006").save()
    assert record.branch_id == 1
    record.delete()
    assert record.branch_id == -1
    # now recreate the same record from ontology_id with a different description
    # there's a unique constraint on ontology_id, so this should recover the trashed record
    record = bt.Ethnicity.from_source(ontology_id="HANCESTRO:0006")
    record.description = "new description"
    record.save()
    # now this record is recovered from the trash with the new description
    assert record.branch_id == 1
    assert record.description == "new description"
    bt.Ethnicity.objects.filter().delete()
