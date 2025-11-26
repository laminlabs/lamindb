import lamindb as ln


def test_view_parents():
    label1 = ln.Record(name="label1")
    label2 = ln.Record(name="label2")
    label1.save()
    label2.save()
    label1.parents.add(label2)
    label1.view_parents(ln.Record.name, distance=1)
    label1.delete(permanent=True)
    label2.delete(permanent=True)


def test_query_parents_children():
    label1 = ln.Record(name="label1").save()
    label2 = ln.Record(name="label2").save()
    label3 = ln.Record(name="label3").save()
    label1.children.add(label2)
    label2.children.add(label3)
    parents = label3.query_parents()
    assert len(parents) == 2
    assert label1 in parents and label2 in parents
    children = label1.query_children()
    assert len(children) == 2
    assert label2 in children and label3 in children
    label1.delete(permanent=True)
    label2.delete(permanent=True)
    label3.delete(permanent=True)


def test_view_lineage_circular():
    import pandas as pd

    transform = ln.Transform(key="test").save()
    run = ln.Run(transform=transform).save()
    artifact = ln.Artifact.from_dataframe(
        pd.DataFrame({"a": [1, 2, 3]}), description="test artifact", run=run
    ).save()
    run.input_artifacts.add(artifact)
    artifact.view_lineage()
    artifact.delete(permanent=True)
    transform.delete(permanent=True)
