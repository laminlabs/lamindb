import lamindb as ln
from lamindb.models.has_parents import _add_emoji


def test_view_parents():
    label1 = ln.ULabel(name="label1")
    label2 = ln.ULabel(name="label2")
    label1.save()
    label2.save()
    label1.parents.add(label2)
    label1.view_parents(ln.ULabel.name, distance=1)
    label1.delete()
    label2.delete()


def test_query_parents_children():
    label1 = ln.ULabel(name="label1").save()
    label2 = ln.ULabel(name="label2").save()
    label3 = ln.ULabel(name="label3").save()
    label1.children.add(label2)
    label2.children.add(label3)
    parents = label3.query_parents()
    assert len(parents) == 2
    assert label1 in parents and label2 in parents
    children = label1.query_children()
    assert len(children) == 2
    assert label2 in children and label3 in children
    label1.delete()
    label2.delete()
    label3.delete()


def test_add_emoji():
    transform = ln.Transform(key="test-12345", type="upload")
    assert _add_emoji(transform, label="transform") == "🖥️ transform"
    transform.save()
    run = ln.Run(transform=transform)
    assert _add_emoji(run, label="run") == "🖥️ run"
    transform.delete()


def test_view_lineage_circular():
    import pandas as pd

    transform = ln.Transform(key="test").save()
    run = ln.Run(transform=transform).save()
    artifact = ln.Artifact.from_df(
        pd.DataFrame({"a": [1, 2, 3]}), description="test artifact", run=run
    ).save()
    run.input_artifacts.add(artifact)
    artifact.view_lineage()
    artifact.delete(permanent=True)
    run.delete()
    transform.delete()
