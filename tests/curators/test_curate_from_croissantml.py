import lamindb as ln


def test_curate_from_croissantml():
    path = ln.examples.croissantml.mini_immuno()
    artifact = ln.integrations.curate_from_croissantml(path)
    assert (
        artifact.description
        == "Mini immuno dataset (mini_immuno.anndata.zarr) - A few samples from the immunology dataset"
    )
    assert artifact.version == "1.0"
    assert (
        artifact.ulabels.get(name="https://creativecommons.org/licenses/by/4.0/")
        is not None
    )
    assert artifact.projects.get(name="Mini Immuno Project") is not None
