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

    record = bt.Ethnicity.from_source(ontology_id="HANCESTRO:0005").save()
    parent = bt.Ethnicity.get(ontology_id="HANCESTRO:0004")
    assert parent in record.parents.list()
    parent.delete()

    bt.Ethnicity.import_from_source()
    parent = bt.Ethnicity.get(ontology_id="HANCESTRO:0004")
    assert parent in record.parents.list()
    record = bt.Ethnicity.get("7RNCY3yC")
    assert record.parents.all().one().name == "South East Asian"
    # the source.in_db should be set to True since we imported all the records
    assert record.source.in_db is True


def test_add_ontology_from_values():
    import bionty as bt

    bt.Ethnicity.filter().delete()
    ln.save(
        bt.Ethnicity.from_values(
            [
                "HANCESTRO:0597",
                "HANCESTRO:0006",
            ],
            field=bt.Ethnicity.ontology_id,
        )
    )
    record = bt.Ethnicity.get("7RNCY3yC")
    assert record.parents.all().one().name == "South East Asian"
    # the source.in_db should be set back to False since we deleted all records
    assert record.source.in_db is False


def test_view_lineage_circular():
    import pandas as pd

    transform = ln.Transform(name="test").save()
    run = ln.Run(transform=transform).save()
    artifact = ln.Artifact.from_df(
        pd.DataFrame({"a": [1, 2, 3]}), description="test artifact", run=run
    ).save()
    run.input_artifacts.add(artifact)
    artifact.view_lineage()
    artifact.delete(permanent=True)
    run.delete()
    transform.delete()
