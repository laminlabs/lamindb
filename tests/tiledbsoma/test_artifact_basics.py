import lamindb as ln
import pytest
from lamindb.models.artifact import data_is_soma_experiment


def test_create_from_soma_experiment(soma_experiment_file, adata_file):
    with pytest.raises(ValueError) as error:
        ln.Artifact.from_tiledbsoma(adata_file, description="test1")
    assert (
        "data has to be a SOMA Experiment object or a path to SOMA Experiment store."
        in error.exconly()
    )

    af = ln.Artifact.from_tiledbsoma(soma_experiment_file, description="test1")
    assert af.description == "test1"
    assert af.key is None
    assert af.otype == "tiledbsoma"
    assert af.n_observations == 3


def test_data_is_soma_experiment_paths():
    assert data_is_soma_experiment("something.tiledbsoma")


def test_data_is_soma_experiment(soma_experiment_file):
    import tiledbsoma

    with tiledbsoma.Experiment.open(soma_experiment_file) as store:
        assert data_is_soma_experiment(store)
