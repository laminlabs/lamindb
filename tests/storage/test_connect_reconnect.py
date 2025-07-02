import lamindb as ln
import pytest


def test_connect_reconnect():
    # testuser2 needs write access lamin-site-assets because of a fluke
    # in the legacy collaborator management, it seems
    assert ln.setup.settings.user.handle == "testuser2"
    ln.connect("lamindb-unit-tests-storage")  # this is not changing anything
    count1 = ln.Artifact.filter().count()
    # a public instance that does not have bionty configured
    ln.connect("laminlabs/lamin-site-assets")
    count2 = ln.Artifact.filter().count()
    assert count1 != count2
    with pytest.raises(ln.setup.errors.ModuleWasntConfigured):
        import bionty as bt
    ln.connect("lamindb-unit-tests-storage")
    import bionty as bt

    count3 = bt.Gene.filter().count()
    assert count2 != count3
