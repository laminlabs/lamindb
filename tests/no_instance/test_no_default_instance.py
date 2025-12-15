import lamindb as ln
import pandas as pd
import pytest


def test_instance_not_connected():
    assert ln.setup.settings.instance.slug == "none/none"

    with pytest.raises(ln.setup.errors.CurrentInstanceNotConfigured):
        ln.Artifact.filter()


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
