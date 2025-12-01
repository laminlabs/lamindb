import os

import lamindb as ln
import pytest


@pytest.mark.skipif(
    os.getenv("LAMINDB_TEST_DB_VENDOR") == "sqlite", reason="Postgres-only"
)
@pytest.mark.parametrize(
    "model_class,model_name",
    [
        (ln.Record, "Record"),
        (ln.ULabel, "ULabel"),
    ],
    ids=["Record", "ULabel"],
)
def test_rename_and_reparent_recordtype(model_class, model_name):
    """Test renaming a record type and changing its parent."""

    # Test simple rename first
    experiment = model_class(name="Experiment", is_type=True).save()
    feature = ln.Feature(name="experiment", dtype=experiment).save()
    experiment.name = "ExperimentRenamed"
    experiment.save()
    feature.refresh_from_db()
    assert feature.dtype == f"cat[{model_name}[ExperimentRenamed]]"

    # Now add a parent (move from root to under a parent)
    parent_type = model_class(name="ParentType", is_type=True).save()
    experiment.type = parent_type
    experiment.save()
    feature.refresh_from_db()
    assert feature.dtype == f"cat[{model_name}[ParentType[ExperimentRenamed]]]"

    # Change to a different parent
    other_parent = model_class(name="OtherParent", is_type=True).save()
    experiment.type = other_parent
    experiment.save()
    feature.refresh_from_db()
    assert feature.dtype == f"cat[{model_name}[OtherParent[ExperimentRenamed]]]"

    # Create a record under the previous parent that has the same name with a feature
    experiment2 = model_class(
        name="ExperimentRenamed", is_type=True, type=parent_type
    ).save()
    feature2 = ln.Feature(name="experiment2", dtype=experiment2).save()
    assert feature2.dtype == f"cat[{model_name}[ParentType[ExperimentRenamed]]]"

    # Test rename the new record type
    experiment2.name = "Experiment"
    experiment2.save()
    feature2.refresh_from_db()
    assert feature2.dtype == f"cat[{model_name}[ParentType[Experiment]]]"
    # this did not mutate the other feature that has the same name
    assert feature.dtype == f"cat[{model_name}[OtherParent[ExperimentRenamed]]]"

    # Remove parent (move back to root)
    experiment.type = None
    experiment.save()
    feature.refresh_from_db()
    assert feature.dtype == f"cat[{model_name}[ExperimentRenamed]]"

    experiment.delete(permanent=True)
    feature.delete(permanent=True)
    experiment2.delete(permanent=True)
    feature2.delete(permanent=True)
    parent_type.delete(permanent=True)
    other_parent.delete(permanent=True)
