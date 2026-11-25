import lamindb as ln
import pytest

FRAME_UID = "gL3TbX2qZQmCwTAU"


@pytest.fixture
def benchmarks_instance():
    """Destination database for the transfer. Outside the timed rounds."""
    ln.connect("laminlabs/lamindb-benchmarks")
    yield


def _delete_transferred() -> None:
    try:
        frame = ln.Record.get(FRAME_UID)
    except ln.Record.DoesNotExist:
        return
    children = frame.records.all()
    if len(children):
        children.delete(permanent=True)
    frame.delete(permanent=True)


def profile_transfer_one_record_frame_with_one_record(benchmark, benchmarks_instance):
    """Time syncing one RNA-seq record frame and one child record."""

    def save():
        ln.core.sync(
            registry=ln.Record,
            uid=FRAME_UID,
            source_db="laminlabs/lamindata",
            depth=1,
            transfer="annotations",
            limit=1,
        )

    def teardown():
        frame = ln.Record.get(FRAME_UID)
        children = frame.records.all()
        assert len(children) == 1, f"Expected 1 child, got {len(children)}"
        children.delete(permanent=True)
        frame.delete(permanent=True)

    benchmark.pedantic(
        save,
        setup=_delete_transferred,
        teardown=teardown,
        rounds=3,
        warmup_rounds=1,
    )
