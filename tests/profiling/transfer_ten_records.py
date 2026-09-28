import lamindb as ln

db = ln.DB("laminlabs/lamindata")
# transfer the RNA-seq frame from lamindata to lamindb-benchmarks
ln.models.sync_objects_from_database(
    registry="record",
    uids="gL3TbX2qZQmCwTAU",
    source="laminlabs/lamindata",
    depth=1,
)


def cleanup():
    transferred_frame = ln.Record.get("gL3TbX2qZQmCwTAU")
    all_children = transferred_frame.records.all()
    assert len(all_children) == 10, f"Expected 10 children, got {len(all_children)}"
    all_children.delete(permanent=True)
    transferred_frame.delete(permanent=True)
