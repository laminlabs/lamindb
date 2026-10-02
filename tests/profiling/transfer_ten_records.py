import lamindb as ln

# transfer the RNA-seq frame from lamindata to lamindb-benchmarks
ln.core.sync(
    registry="record",
    uid="gL3TbX2qZQmCwTAU",
    source_db="laminlabs/lamindata",
    depth=1,
    transfer="annotations",
)


def cleanup():
    transferred_frame = ln.Record.get("gL3TbX2qZQmCwTAU")
    all_children = transferred_frame.records.all()
    assert len(all_children) == 10, f"Expected 10 children, got {len(all_children)}"
    all_children.delete(permanent=True)
    transferred_frame.delete(permanent=True)
