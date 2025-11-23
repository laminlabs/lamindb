import lamindb as ln


def test_simple_rename_no_parent():
    """Test renaming a record type with no parent (root level)."""

    cell_type = ln.Record(name="CellType", is_type=True).save()

    feature = ln.Feature(name="cell_annotation", dtype=cell_type).save()

    cell_type.name = "CellTypeRenamed"
    cell_type.save()

    feature.refresh_from_db()

    assert feature.dtype == "cat[Record[CellTypeRenamed]]"
