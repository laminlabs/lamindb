import uuid

import lamindb as ln
import pytest

ARTIFACT_KEY = "profile/mini_immuno.parquet"


@pytest.fixture
def mini_immuno(tmp_path):
    """Local instance, features, and the 3-row mini-immuno dataframe.

    Setup and teardown stay outside the timed `benchmark.pedantic` rounds.
    """
    instance_name = f"curate-profile-{uuid.uuid4().hex[:8]}"
    ln.setup.init(storage=tmp_path / "storage", name=instance_name, modules="bionty")
    try:
        ln.examples.datasets.mini_immuno.define_features_labels()
        schema = ln.examples.datasets.mini_immuno.define_mini_immuno_schema_flexible()
        df = ln.examples.datasets.mini_immuno.get_dataset1(otype="DataFrame")
        yield df, schema
    finally:
        ln.setup.delete(instance_name, force=True)


def profile_curate_simple_dataframe(benchmark, mini_immuno):
    """Time Artifact.from_dataframe().save() on the mini-immuno dataframe."""
    df, schema = mini_immuno

    def save():
        ln.Artifact.from_dataframe(df, key=ARTIFACT_KEY, schema=schema).save()

    def teardown():
        ln.Artifact.get(key=ARTIFACT_KEY).delete(permanent=True)

    benchmark.pedantic(save, teardown=teardown, rounds=5, iterations=1, warmup_rounds=1)
