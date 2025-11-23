import lamindb as ln


def test_rename_and_reparent_recordtype():
    """Test renaming a record type and changing its parent."""

    # Test simple rename first
    experiment = ln.Record(name="Experiment", is_type=True).save()
    feature = ln.Feature(name="cell_annotation", dtype=experiment).save()
    experiment.name = "ExperimentRenamed"
    experiment.save()
    feature.refresh_from_db()
    assert feature.dtype == "cat[Record[ExperimentRenamed]]"

    # Now add a parent (move from root to under a parent)
    parent_type = ln.Record(name="ParentType", is_type=True).save()
    experiment.type = parent_type
    experiment.save()
    feature.refresh_from_db()
    assert feature.dtype == "cat[Record[ParentType[ExperimentRenamed]]]"

    # Change to a different parent
    other_parent = ln.Record(name="OtherParent", is_type=True).save()
    experiment.type = other_parent
    experiment.save()
    feature.refresh_from_db()
    assert feature.dtype == "cat[Record[OtherParent[ExperimentRenamed]]]"

    # Remove parent (move back to root)
    experiment.type = None
    experiment.save()
    feature.refresh_from_db()
    assert feature.dtype == "cat[Record[ExperimentRenamed]]"
