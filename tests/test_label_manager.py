import lnschema_bionty as lb
import pytest

import lamindb as ln


@pytest.fixture
def get_test_files():
    with open("./default_storage/test-inherit1", "w") as f:
        f.write("file1")
    with open("./default_storage/test-inherit2", "w") as f:
        f.write("file2")
    file1 = ln.File("./default_storage/test-inherit1")
    file1.save()
    file2 = ln.File("./default_storage/test-inherit2")
    file2.save()
    yield file1, file2
    file1.delete(permanent=True, storage=True)
    file2.delete(permanent=True, storage=True)


# also see test_feature_manager!
def test_add_from(get_test_files):
    file1, file2 = get_test_files
    label_names = [f"Project {i}" for i in range(3)]
    labels = [ln.ULabel(name=label_name) for label_name in label_names]
    ln.save(labels)

    cell_line_names = [f"Cell line {i}" for i in range(3)]
    cell_lines = [lb.CellLine(name=name) for name in cell_line_names]
    ln.save(cell_lines)

    file2.ulabels.add(*labels)
    # here test add without passing a feature
    file2.labels.add(cell_lines)
    assert file2.cell_lines.count() == len(cell_lines)

    assert file1.ulabels.exists() is False
    file1.labels.add_from(file2)
    assert file1.ulabels.count() == file2.ulabels.count()
    assert file1.cell_lines.count() == file2.cell_lines.count()

    for label in labels:
        label.delete()
    for cell_line in cell_lines:
        cell_line.delete()
