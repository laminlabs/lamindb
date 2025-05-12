import concurrent.futures

import lamindb as ln
import pandas as pd
import pytest


@ln.tracked()
def process_chunk(chunk_id: int) -> str:
    # Create a simple DataFrame
    df = pd.DataFrame(
        {"id": range(chunk_id * 10, (chunk_id + 1) * 10), "value": range(10)}
    )

    # Save it as an artifact
    key = f"chunk_{chunk_id}.parquet"
    artifact = ln.Artifact.from_df(df, key=key).save()
    return artifact.key


def test_tracked_parallel():
    param_type = ln.Feature(name="Script[test_tracked.py]", is_type=True).save()
    ln.Feature(name="chunk_id", dtype="int", type=param_type).save()

    with pytest.raises(RuntimeError) as err:
        process_chunk(4)
    assert (
        err.exconly()
        == "RuntimeError: Please track the global run context before using @ln.tracked(): ln.track()"
    )

    # Ensure tracking is on
    ln.track()

    # Number of parallel executions
    n_parallel = 3

    # Use ThreadPoolExecutor for parallel execution
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_parallel) as executor:
        # Submit all tasks
        futures = [executor.submit(process_chunk, i) for i in range(n_parallel)]
        # Get results as they complete
        chunk_keys = [
            future.result() for future in concurrent.futures.as_completed(futures)
        ]

    # Verify results
    # Each execution should have created its own artifact with unique run
    print(f"Created artifacts with keys: {chunk_keys}")
    artifacts = [ln.Artifact.get(key=key) for key in chunk_keys]

    # Check that we got the expected number of artifacts
    assert len(artifacts) == n_parallel

    # Verify each artifact has its own unique run
    runs = [artifact.run for artifact in artifacts]
    run_ids = [run.id for run in runs]
    print(f"Run IDs: {run_ids}")
    assert len(set(run_ids)) == n_parallel  # all runs should be unique

    # Verify each run has the correct start and finish times
    for run in runs:
        print(f"Run details: {run}")
        assert run.started_at is not None
        assert run.finished_at is not None
        assert run.started_at < run.finished_at

    # Clean up test artifacts
    for artifact in artifacts:
        artifact.delete(permanent=True)

    ln.context._uid = None
    ln.context._run = None
    ln.context._transform = None
    ln.context._path = None


if __name__ == "__main__":
    test_tracked_parallel()
