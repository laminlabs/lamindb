import lnschema_bionty as lb
import pytest

import lamindb as ln
from lamindb._parents import _add_emoji


def test_view_parents():
    label1 = ln.Label(name="label1")
    label2 = ln.Label(name="label2")
    label1.save()
    label2.save()
    label1.parents.add(label2)
    label1.view_parents(ln.Label.name, distance=1)
    label1.delete()
    label2.delete()

    with pytest.raises(NotImplementedError):
        lb.Species(name="mouse").view_parents()


def add_emoji():
    assert _add_emoji(ln.Transform(type="app"), label="transform") == "🖥️ transform"
    assert (
        _add_emoji(ln.Run(transform=ln.Transform(type="app")), label="run") == "🖥️ run"
    )
