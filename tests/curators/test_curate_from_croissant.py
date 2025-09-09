import shutil

import lamindb as ln
import pytest


@pytest.mark.parametrize("filepath_prefix", [None, "testdb/"])
def test_curate_artifact_from_croissant(filepath_prefix: str | None):
    croissant_path, dataset1_path = ln.examples.croissant.mini_immuno(
        n_files=1, filepath_prefix=filepath_prefix
    )
    artifact1 = ln.integrations.curate_from_croissant(croissant_path)
    croissant_path.unlink()
    shutil.rmtree(dataset1_path)
    assert (
        artifact1.description
        == "Mini immuno dataset - A few samples from the immunology dataset"
    )
    assert artifact1.key == "mini_immuno.anndata.zarr"
    assert artifact1.version == "1.0"
    assert (
        artifact1._key_is_virtual
        if filepath_prefix is None
        else not artifact1._key_is_virtual
    )
    license_label = artifact1.ulabels.get(
        name="https://creativecommons.org/licenses/by/4.0/"
    )
    project_label = artifact1.projects.get(name="Mini Immuno Project")

    artifact1.delete(permanent=True, storage=True)  # because of real storage key
    project_label.delete(permanent=True)
    license_label.delete(permanent=True)


def test_curate_collection_from_croissant():
    croissant_path, dataset1_path, dataset2_path = ln.examples.croissant.mini_immuno(
        n_files=2
    )
    collection = ln.integrations.curate_from_croissant(croissant_path)
    croissant_path.unlink()
    shutil.rmtree(dataset1_path)
    dataset2_path.unlink()
    artifact1 = collection.artifacts.get(key="mini_immuno.anndata.zarr")
    artifact2 = collection.artifacts.get(key="mini.csv")
    license_label = collection.ulabels.get(
        name="https://creativecommons.org/licenses/by/4.0/"
    )
    project_label = collection.projects.get(name="Mini Immuno Project")

    collection.delete(permanent=True)
    artifact1.delete(permanent=True)
    artifact2.delete(permanent=True)
    project_label.delete(permanent=True)
    license_label.delete(permanent=True)
