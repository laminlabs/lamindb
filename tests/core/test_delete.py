import lamindb as ln
import pytest


@pytest.mark.parametrize("permanent", [True, False])
@pytest.mark.parametrize("size", [3, 10000])
def test_delete_qs(permanent, size):
    """Test deletion behavior for small (<10000) and large (>=10000) querysets.

    Small querysets delete individually, large ones trigger bulk delete."""
    prefix = f"testlabel_{size}"
    ln.settings.creation.search_names = False
    labels = [ln.ULabel(name=f"{prefix}{i}") for i in range(size)]
    ln.settings.creation.search_names = True
    ln.save(labels)
    ln.ULabel.filter(name__startswith=prefix).delete(permanent=permanent)
    assert ln.ULabel.filter(name__startswith=prefix, branch_id=-1).count() == (
        0 if permanent else size
    )
    assert ln.ULabel.filter(name__startswith=prefix).count() == 0
