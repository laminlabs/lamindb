import lamindb as ln
import pytest


def test_connect_reconnect():
    ln.connect("laminlabs/lamindata")
    count1 = ln.Artifact.filter().count()
    ln.connect("laminlabs/lamin-site-assets")
    count2 = ln.Artifact.filter().count()
    assert count1 != count2
    with pytest.raises(ln.setup.errors.ModuleWasntConfigured):
        import bionty as bt
    ln.connect("laminlabs/lamindata")
    import bionty as bt

    count3 = bt.Gene.filter().count()
    assert count2 != count3
