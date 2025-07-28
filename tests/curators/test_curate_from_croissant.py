import shutil

import lamindb as ln


def test_curate_from_croissant():
    croissant_path, dataset1_path = ln.examples.croissant.mini_immuno(n_files=1)
    artifact1 = ln.integrations.curate_from_croissant(croissant_path)
    croissant_path.unlink()
    shutil.rmtree(dataset1_path)
    assert (
        artifact1.description
        == "Mini immuno dataset (mini_immuno.anndata.zarr) - A few samples from the immunology dataset"
    )
    assert artifact1.version == "1.0"
    license_label = artifact1.ulabels.get(
        name="https://creativecommons.org/licenses/by/4.0/"
    )
    project_label = artifact1.projects.get(name="Mini Immuno Project")
    artifact1.delete(permanent=True)
    project_label.delete()
    license_label.delete()
