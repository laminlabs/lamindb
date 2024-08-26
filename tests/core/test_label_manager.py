import bionty as bt
import lamindb as ln
import pytest


@pytest.fixture
def get_test_artifacts():
    with open("./default_storage_unit_core/test-inherit1", "w") as f:
        f.write("artifact1")
    with open("./default_storage_unit_core/test-inherit2", "w") as f:
        f.write("artifact2")
    artifact1 = ln.Artifact("./default_storage_unit_core/test-inherit1")
    artifact1.save()
    artifact2 = ln.Artifact("./default_storage_unit_core/test-inherit2")
    artifact2.save()
    yield artifact1, artifact2
    artifact1.delete(permanent=True, storage=True)
    artifact2.delete(permanent=True, storage=True)


# also see test_feature_manager!
def test_add_from(get_test_artifacts):
    artifact1, artifact2 = get_test_artifacts
    label_names = [f"Project {i}" for i in range(3)]
    ulabels = [ln.ULabel(name=label_name) for label_name in label_names]
    ln.save(ulabels)

    cell_line_names = [f"Cell line {i}" for i in range(3)]
    cell_lines = [bt.CellLine(name=name) for name in cell_line_names]
    ln.save(cell_lines)

    # pass a list of length 0
    artifact2.labels.add([])
    # now actually pass the labels
    artifact2.labels.add(ulabels)
    # here test add without passing a feature
    artifact2.labels.add(cell_lines)
    assert artifact2.cell_lines.count() == len(cell_lines)

    assert artifact1.ulabels.exists() is False
    artifact1.labels.add_from(artifact2)
    assert artifact1.ulabels.count() == artifact2.ulabels.count()
    assert artifact1.cell_lines.count() == artifact2.cell_lines.count()

    artifact2.cell_lines.remove(*cell_lines)
    artifact1.cell_lines.remove(*cell_lines)
    artifact2.ulabels.remove(*ulabels)
    artifact1.ulabels.remove(*ulabels)

    for ulabel in ulabels:
        ulabel.delete()
    for cell_line in cell_lines:
        cell_line.delete()
