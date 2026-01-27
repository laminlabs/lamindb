from pathlib import Path
from typing import Iterable

import lamindb as ln
import pandas as pd
import pytest


@ln.flow(global_run="clear")
def process_chunk(
    chunk_id: int, artifact_param: ln.Artifact, records_params: Iterable[ln.Record]
) -> str:
    # Create a simple DataFrame
    df = pd.DataFrame(
        {"id": range(chunk_id * 10, (chunk_id + 1) * 10), "value": range(10)}
    )
    env_file = Path("file_with_same_hash.txt")
    env_file.write_text("1")
    ln.Artifact(env_file, description="file_with_same_hash").save()
    # Save it as an artifact
    key = f"chunk_{chunk_id}.parquet"
    artifact = ln.Artifact.from_dataframe(df, key=key).save()
    assert ln.context.run is not None
    return artifact.key


def test_flow():
    param_artifact = ln.Artifact(".gitignore", key="param_artifact").save()
    ln.Record(name="record1").save(), ln.Record(name="record2").save()
    records_params = ln.Record.filter(name__startswith="record")

    assert ln.context.run is None
    artifact_key = process_chunk(1, param_artifact, records_params)
    assert ln.context.run is None

    # Verify the artifacts and runs
    artifacts = [ln.Artifact.get(key=key) for key in [artifact_key]]
    same_hash_artifacts = ln.Artifact.filter(description="file_with_same_hash")

    runs = [artifact.run for artifact in artifacts]

    # Verify each run has the correct start and finish times
    for run in runs:
        print(f"Run details: {run}")
        assert run.started_at is not None
        assert run.finished_at is not None
        assert run.started_at < run.finished_at
        assert run.status == "completed"
        assert isinstance(run.params["chunk_id"], int)
        assert run.params["artifact_param"].startswith(
            f"Artifact[{param_artifact.uid}]"
        )
        assert run.params["records_params"] == [
            f"Record[{record.uid}]" for record in records_params
        ]

    # test error behavior
    with pytest.raises(RuntimeError) as error:
        ln.context._run = run
        process_chunk(1, param_artifact, records_params)
        ln.context._run = None
    assert str(error.exconly()).startswith(
        "RuntimeError: Please clear the global run context before using @ln.flow(): no `ln.track()` or `@ln.flow(global_run='clear')`"
    )

    # Clean up test artifacts
    runs = []
    for artifact in artifacts:
        runs.append(artifact.run)
        artifact.delete(permanent=True)
    param_artifact.delete(permanent=True)
    same_hash_artifacts[0].delete(permanent=True)
    Path("file_with_same_hash.txt").unlink()
    for run in runs:
        run.delete(permanent=True)
