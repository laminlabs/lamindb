import time
from pathlib import Path
from typing import Iterable

import lamindb as ln
import pandas as pd
import pytest
from lamindb.errors import InvalidArgument


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
        "RuntimeError: Please use @ln.step() or clear the global run context before using @ln.flow(): no `ln.track()` or `@ln.flow(global_run='clear')`"
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
    ln.context._run = None


def test_flow_track_arg_aliases_implicit():
    unique = time.time_ns()
    missing_project = f"missing-flow-project-{unique}"

    @ln.flow(global_run="clear")
    def flow_with_implicit_project_alias(project: str) -> None:
        pass

    with pytest.raises(InvalidArgument) as error:
        flow_with_implicit_project_alias(project=missing_project)
    assert error.exconly().startswith(
        f"lamindb.errors.InvalidArgument: Project '{missing_project}' not found"
    )


def test_flow_track_arg_aliases_false():
    unique = time.time_ns()
    missing_project = f"missing-flow-project-{unique}"

    @ln.flow(global_run="clear", track_arg_aliases=False)
    def flow_without_project_alias(project: str) -> str:
        assert ln.context.run is not None
        return ln.context.run.uid

    run = None
    try:
        run_uid = flow_without_project_alias(project=missing_project)
        run = ln.Run.get(uid=run_uid)
        assert run.params["project"] == missing_project
    finally:
        ln.context._run = None
        if run is not None:
            run.delete(permanent=True)
            run.transform.delete(permanent=True)


def test_flow_exception_uploads_run_report():
    run_uid: str | None = None
    run = None
    transform = None
    report = None
    try:

        @ln.flow(global_run="clear")
        def failing_flow() -> None:
            nonlocal run_uid
            assert ln.context.run is not None
            run_uid = ln.context.run.uid
            print("before flow error")
            raise ValueError("flow failed")

        with pytest.raises(ValueError, match="flow failed"):
            failing_flow()

        assert ln.context.run is None
        assert run_uid is not None
        run = ln.Run.get(uid=run_uid)
        transform = run.transform
        report = run.report
        assert run.status == "errored"
        assert report is not None
        report_text = report.cache().read_text()
        assert "before flow error" in report_text
    finally:
        ln.context._run = None
        if run is not None:
            run.delete(permanent=True)
        if report is not None:
            report.delete(permanent=True)
        if transform is not None:
            transform.delete(permanent=True)


def test_flow_annotation_mismatch_skips_param():
    run = None
    try:

        @ln.flow(global_run="clear")
        def typed_flow(count: int, label: str) -> str:
            assert ln.context.run is not None
            return ln.context.run.uid

        run_uid = typed_flow("not-an-int", "ok")
        run = ln.Run.get(uid=run_uid)
        assert run.params == {"label": "ok"}
    finally:
        ln.context._run = None
        if run is not None:
            run.delete(permanent=True)


def test_flow_annotation_match_keeps_param():
    run = None
    try:

        @ln.flow(global_run="clear")
        def typed_flow(count: int, label: str) -> str:
            assert ln.context.run is not None
            return ln.context.run.uid

        run_uid = typed_flow(3, "ok")
        run = ln.Run.get(uid=run_uid)
        assert run.params == {"count": 3, "label": "ok"}
    finally:
        ln.context._run = None
        if run is not None:
            run.delete(permanent=True)


# see test_record_features_add_remove_values in test_record_basics.py
# for a similar broad dtype matrix on Record features
def test_flow_annotation_runtime_validation_type_matrix():
    run_valid = None
    run_invalid = None
    record = None
    record2 = None
    try:
        unique = time.time_ns()
        record = ln.Record(name=f"track-param-record-{unique}").save()
        record2 = ln.Record(name=f"track-param-record2-{unique}").save()

        @ln.flow(global_run="clear")
        def typed_flow(
            bool_value: bool,
            str_value: str,
            list_str_value: list[str],
            int_value: int,
            list_int_value: list[int],
            float_value: float,
            list_float_value: list[float],
            dict_value: dict,
            record_value: ln.Record,
            records_value: list[ln.Record],
        ) -> str:
            assert ln.context.run is not None
            return ln.context.run.uid

        run_uid = typed_flow(
            True,
            "alpha",
            ["a", "b"],
            42,
            [1, 2, 3],
            3.14,
            [2.71, 3.14],
            {"key": "value"},
            record,
            [record, record2],
        )
        run_valid = ln.Run.get(uid=run_uid)
        assert run_valid.params == {
            "bool_value": True,
            "str_value": "alpha",
            "list_str_value": ["a", "b"],
            "int_value": 42,
            "list_int_value": [1, 2, 3],
            "float_value": 3.14,
            "list_float_value": [2.71, 3.14],
            "dict_value": {"key": "value"},
            "record_value": f"Record[{record.uid}]",
            "records_value": [f"Record[{record.uid}]", f"Record[{record2.uid}]"],
        }

        run_uid = typed_flow(
            "not-a-bool",
            "alpha",
            ["a", "b"],
            "not-an-int",
            [1, 2, 3],
            3.14,
            [2.71, 3.14],
            {"key": "value"},
            "not-a-record",
            ["not-a-record"],
        )
        run_invalid = ln.Run.get(uid=run_uid)
        assert run_invalid.params == {
            "str_value": "alpha",
            "list_str_value": ["a", "b"],
            "list_int_value": [1, 2, 3],
            "float_value": 3.14,
            "list_float_value": [2.71, 3.14],
            "dict_value": {"key": "value"},
        }
    finally:
        ln.context._run = None
        if run_valid is not None:
            run_valid.delete(permanent=True)
        if run_invalid is not None:
            run_invalid.delete(permanent=True)
        if record is not None:
            record.delete(permanent=True)
        if record2 is not None:
            record2.delete(permanent=True)
