import lamindb as ln
import pytest


@pytest.mark.parametrize("permanent", [True, False])
def test_delete_qs(permanent):
    """Test deletion behavior for small (1) and large (>=2) querysets.

    Small querysets delete individually, large ones trigger bulk delete."""
    ln.settings.creation.search_names = False
    labels = [ln.ULabel(name=f"label_{i}") for i in range(3)]
    ln.settings.creation.search_names = True
    ln.save(labels)
    ln.ULabel.filter(name__startswith="label_").delete(permanent=permanent)
    assert ln.ULabel.filter(name__startswith="label_", branch_id=-1).count() == (
        0 if permanent else 3
    )
    assert ln.ULabel.filter(name__startswith="label_").count() == 0
