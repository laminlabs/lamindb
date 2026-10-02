import bionty as bt
import lamindb as ln
import pytest


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


def test_transfer_run_label_uses_description():
    from lamindb.models.has_parents import get_record_label

    transform = ln.Transform(
        key="__lamindb_transfer__/4XIuR0tvaiXM",
        description="Transfer from `laminlabs/lamindata`",
        kind="function",
        uid="4XIuR0tvaiXM0000",
    ).save()
    run = ln.Run(transform=transform).save()
    label = get_record_label(run)
    assert "Transfer from `laminlabs/lamindata`" in label
    assert "__lamindb_transfer__" not in label

    script = ln.Transform(key="my-script.py").save()
    script_run = ln.Run(transform=script).save()
    assert "my-script.py" in get_record_label(script_run)

    run.delete(permanent=True)
    script_run.delete(permanent=True)
    transform.delete(permanent=True)
    script.delete(permanent=True)


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


def test_view_parents_connected_instance():
    ct = bt.CellType.connect("laminlabs/cellxgene").first()

    if ct and hasattr(ct, "parents"):
        ct.view_parents(distance=2, with_children=True)


def test_query_relatives_connected_instance():
    ct = bt.CellType.connect("laminlabs/cellxgene").filter(name="T cell").first()

    if ct:
        parents = ct.query_parents()
        assert parents.db == "laminlabs/cellxgene"

        children = ct.query_children()
        assert children.db == "laminlabs/cellxgene"


def test_view_lineage_connected_instance():
    af = ln.Artifact.connect("laminlabs/cellxgene").first()

    if af and af.run:
        af.view_lineage()


@pytest.mark.parametrize("terminal_ipython", [False, True])
def test_view_digraph_keeps_rendered_files_out_of_working_directory(
    tmp_path, monkeypatch, terminal_ipython
):
    from pathlib import Path

    import graphviz
    from lamindb.models import has_parents

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(has_parents, "is_run_from_ipython", terminal_ipython)
    if terminal_ipython:
        import IPython

        terminal_shell = type("TerminalInteractiveShell", (), {})()
        monkeypatch.setattr(IPython, "get_ipython", lambda: terminal_shell)
    viewed = []
    monkeypatch.setattr(
        graphviz.Digraph,
        "_view",
        lambda self, path, **kwargs: viewed.append(path),
    )
    graph = graphviz.Digraph("regression-lineage")
    graph.edge("input", "output")

    rendered = Path(has_parents.view_digraph(graph))

    assert not list(tmp_path.iterdir())
    assert rendered.is_file()
    assert rendered.read_bytes().startswith(b"%PDF")
    assert not rendered.with_suffix("").exists()
    assert viewed == [str(rendered)]


def test_view_digraph_notebook_does_not_create_files(tmp_path, monkeypatch):
    import graphviz
    import IPython
    import IPython.display
    from lamindb.models import has_parents

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(has_parents, "is_run_from_ipython", True)
    notebook_shell = type("ZMQInteractiveShell", (), {})()
    monkeypatch.setattr(IPython, "get_ipython", lambda: notebook_shell)
    displayed = []
    monkeypatch.setattr(
        IPython.display,
        "display",
        lambda bundle, **kwargs: displayed.append(bundle),
    )
    graph = graphviz.Digraph("notebook-lineage")
    graph.edge("input", "output")

    assert has_parents.view_digraph(graph) is None
    assert not list(tmp_path.iterdir())
    assert "image/svg+xml" in displayed[0]
