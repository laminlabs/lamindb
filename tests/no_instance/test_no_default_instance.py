import lamindb as ln
import pandas as pd
import pytest
from lamindb_setup.errors import CurrentInstanceNotConfigured


def test_no_track_run_input_warning_without_global_instance(ccaplog):
    from lamindb.models.artifact import WARNING_NO_INPUT, track_run_input

    assert ln.setup.settings.instance.slug == "none/none"

    track_run_input([])

    assert WARNING_NO_INPUT not in ccaplog.text


def test_no_read_only_warning(ccaplog):
    ln.Artifact.connect("laminlabs/lamindata")
    ln.DB("laminlabs/lamindata")

    assert "connected in read-only mode" not in ccaplog.text


def test_instance_not_connected():
    assert ln.setup.settings.instance.slug == "none/none"

    with pytest.raises(CurrentInstanceNotConfigured):
        ln.Artifact.filter().count()


def test_query_artifacts_lamindata():
    artifacts = ln.Artifact.connect("laminlabs/lamindata")
    n_artifacts = artifacts.count()
    assert n_artifacts > 0
    assert n_artifacts > artifacts.filter().count()


def test_get_artifact_lamindata():
    artifact = ln.Artifact.connect("laminlabs/lamindata").get(
        key="example_datasets/small_dataset1.parquet"
    )
    assert isinstance(artifact.load(), pd.DataFrame)
