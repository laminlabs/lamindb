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

    bt.CellType.from_values(["B cell", "T cell"]).save()
    assert bt.CellType.get(name="B cell").parents.to_dataframe().shape[0] == 1
    bt.CellType.filter().delete(permanent=True)


def test_save_batch_size():
    label_names = [f"Record {i} batch_size" for i in range(3)]
    labels = [ln.Record(name=name) for name in label_names]
    # test bulk creation of new records with batch size
    ln.save(labels, batch_size=2)
    assert ln.Record.filter(name__in=label_names).distinct().count() == 3


def test_bulk_save_lazy_record_features():
    cell_type = ln.Record(name="lazy-cell-type", is_type=True).save()
    ln.Record(name="lazy-b-cell", type=cell_type).save()
    ln.Record(name="lazy-t-cell", type=cell_type).save()
    score_feature = ln.Feature(name="lazy-bulk-score", dtype=float).save()
    cell_feature = ln.Feature(name="lazy-bulk-cell", dtype=cell_type).save()
    schema = ln.Schema([score_feature, cell_feature], name="lazy-bulk-schema").save()
    sheet = ln.Record(name="lazy-sheet", is_type=True, schema=schema).save()

    records = [
        ln.Record(
            name="lazy-sample-1",
            type=sheet,
            features={"lazy-bulk-score": 0.1, "lazy-bulk-cell": "lazy-b-cell"},
        ),
        ln.Record(
            name="lazy-sample-2",
            type=sheet,
            features={"lazy-bulk-score": 0.2, "lazy-bulk-cell": "lazy-t-cell"},
        ),
    ]
    ln.save(records)

    sample_1 = ln.Record.get(name="lazy-sample-1")
    sample_2 = ln.Record.get(name="lazy-sample-2")
    sample_1_values = sample_1.features.get_values()
    sample_2_values = sample_2.features.get_values()
    assert sample_1_values["lazy-bulk-score"] == 0.1
    assert sample_2_values["lazy-bulk-score"] == 0.2
    assert sample_1_values["lazy-bulk-cell"] == "lazy-b-cell"
    assert sample_2_values["lazy-bulk-cell"] == "lazy-t-cell"
    assert not hasattr(records[0], "_features")
    assert not hasattr(records[1], "_features")

    ln.Record.filter(name__in=["lazy-sample-1", "lazy-sample-2"]).delete(permanent=True)
    ln.Record.filter(name="lazy-sheet").delete(permanent=True)
    ln.Record.filter(name__in=["lazy-b-cell", "lazy-t-cell"]).delete(permanent=True)
    ln.Record.filter(name="lazy-cell-type").delete(permanent=True)
    schema.delete(permanent=True)
    score_feature.delete(permanent=True)
    cell_feature.delete(permanent=True)


def test_bulk_save_lazy_record_features_requires_same_schema():
    feature_a = ln.Feature(name="lazy-schema-a", dtype=float).save()
    feature_b = ln.Feature(name="lazy-schema-b", dtype=float).save()
    schema_a = ln.Schema([feature_a], name="lazy-schema-a").save()
    schema_b = ln.Schema([feature_b], name="lazy-schema-b").save()
    type_a = ln.Record(name="lazy-type-a", is_type=True, schema=schema_a).save()
    type_b = ln.Record(name="lazy-type-b", is_type=True, schema=schema_b).save()

    records = [
        ln.Record(name="lazy-mixed-1", type=type_a, features={"lazy-schema-a": 1.0}),
        ln.Record(name="lazy-mixed-2", type=type_b, features={"lazy-schema-b": 2.0}),
    ]
    with pytest.raises(
        ln.errors.ValidationError,
        match="same type schema",
    ):
        ln.save(records)

    ln.Record.filter(name__in=["lazy-mixed-1", "lazy-mixed-2"]).delete(permanent=True)
    ln.Record.filter(name__in=["lazy-type-a", "lazy-type-b"]).delete(permanent=True)
    schema_a.delete(permanent=True)
    schema_b.delete(permanent=True)
    feature_a.delete(permanent=True)
    feature_b.delete(permanent=True)


def test_bulk_save_lazy_record_features_requires_schema():
    unschematized_type = ln.Record(name="lazy-no-schema-type", is_type=True).save()

    records = [
        ln.Record(
            name="lazy-no-schema-1", type=unschematized_type, features={"foo": 1.0}
        )
    ]
    with pytest.raises(
        ln.errors.ValidationError,
        match="same non-null type schema",
    ):
        ln.save(records)

    ln.Record.filter(name="lazy-no-schema-1").delete(permanent=True)
    ln.Record.filter(name="lazy-no-schema-type").delete(permanent=True)


def test_bulk_resave_trashed_records():
    import bionty as bt

    # first create records from public source
    records = bt.Ethnicity.from_values(["asian", "white"]).save()
    assert len(records) == 2
    # parents are also created
    ethnicities = bt.Ethnicity.filter()
    assert ethnicities.count() > 2
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
    assert ethnicities.count() > 3

    # clean up
    ethnicities.delete(permanent=True)
