import shutil

import lamindb as ln


def test_curate_from_croissantml():
    dataset_path, croissantml_path = ln.examples.croissantml.mini_immuno()
    artifact = ln.integrations.curate_from_croissantml(croissantml_path)
    shutil.rmtree(dataset_path)
    assert (
        artifact.description
        == "Mini immuno dataset (mini_immuno.anndata.zarr) - A few samples from the immunology dataset"
    )
    assert artifact.version == "1.0"
    license_label = artifact.ulabels.get(
        name="https://creativecommons.org/licenses/by/4.0/"
    )
    project_label = artifact.projects.get(name="Mini Immuno Project")
    artifact.delete(permanent=True)
    project_label.delete()
    license_label.delete()
