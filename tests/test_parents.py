import lamindb as ln
from lamindb._parents import _add_emoji


def test_view_parents():
    label1 = ln.ULabel(name="label1")
    label2 = ln.ULabel(name="label2")
    label1.save()
    label2.save()
    label1.parents.add(label2)
    label1.view_parents(ln.ULabel.name, distance=1)
    label1.delete()
    label2.delete()


def test_add_emoji():
    record = ln.Transform(type="app")
    assert _add_emoji(record, label="transform") == "🖥️ transform"
    transform = ln.Transform(name="test", type="app")
    transform.save()
    record = ln.Run(transform=transform)
    assert _add_emoji(record, label="run") == "🖥️ run"
    transform.delete()


def test_add_ontology_from_df():
    import bionty as bt

    bt.Ethnicity.filter().delete()
    bt.Ethnicity.save_from_df()
    record = bt.Ethnicity.get("7RNCY3yC")
    assert record.parents.all().one().name == "South East Asian"
