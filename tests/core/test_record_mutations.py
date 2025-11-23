import lamindb as ln


def test_rename_and_reparent_recordtype():
    """Test renaming a record type and changing its parent."""

    # Test simple rename first
    experiment = ln.Record(name="Experiment", is_type=True).save()
    feature = ln.Feature(name="experiment", dtype=experiment).save()
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

    # Create a record under the previous parent that has the same name with a feature
    experiment2 = ln.Record(
        name="ExperimentRenamed", is_type=True, type=parent_type
    ).save()
    feature2 = ln.Feature(name="experiment2", dtype=experiment2).save()
    assert feature2.dtype == "cat[Record[ParentType[ExperimentRenamed]]]"

    # Test rename the new record type
    experiment2.name = "Experiment"
    experiment2.save()
    feature2.refresh_from_db()
    assert feature2.dtype == "cat[Record[ParentType[Experiment]]]"
    # this did not mutate the other feature that has the same name
    assert feature.dtype == "cat[Record[OtherParent[ExperimentRenamed]]]"

    # Remove parent (move back to root)
    experiment.type = None
    experiment.save()
    feature.refresh_from_db()
    assert feature.dtype == "cat[Record[ExperimentRenamed]]"

    experiment.delete(permanent=True)
    feature.delete(permanent=True)
    experiment2.delete(permanent=True)
    feature2.delete(permanent=True)
    parent_type.delete(permanent=True)
    other_parent.delete(permanent=True)
