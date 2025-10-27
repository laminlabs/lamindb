import shutil

import lamindb as ln
import pytest


@pytest.mark.parametrize("filepath_prefix", [None, "test-curators-db/"])
def test_curate_artifact_from_croissant(filepath_prefix: str | None):
    croissant_path, dataset1_path = ln.examples.croissant.mini_immuno(
        n_files=1, filepath_prefix=filepath_prefix
    )
    artifact1 = ln.integrations.curate_from_croissant(croissant_path)
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

    # now mutate the dataset and create a new version
    croissant_path, dataset1_path = ln.examples.croissant.mini_immuno(
        n_files=1, filepath_prefix=filepath_prefix, strip_version=True
    )
    dummy_file_path = dataset1_path / "dummy_file.txt"
    dummy_file_path.write_text("dummy file")

    artifact2 = ln.integrations.curate_from_croissant(croissant_path)
    assert artifact2.description == artifact1.description
    assert artifact2.key == artifact1.key
    assert artifact2.version is None
    assert artifact2.stem_uid == artifact1.stem_uid
    assert artifact2.uid != artifact1.uid
    assert (
        artifact2._key_is_virtual
        if filepath_prefix is None
        else not artifact1._key_is_virtual
    )
    license_label = artifact2.ulabels.get(
        name="https://creativecommons.org/licenses/by/4.0/"
    )
    project_label = artifact2.projects.get(name="Mini Immuno Project")

    shutil.rmtree(dataset1_path)
    croissant_path.unlink()
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
