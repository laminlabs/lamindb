import re

import lamindb as ln
import pandas as pd
from lamindb.examples.fixtures.sheets import (
    populate_nextflow_sheet_with_samples,  # noqa: F401
    populate_sheets_compound_treatment,  # noqa: F401
)


def _assert_describe_feature_columns(
    describe_str: str,
    n_columns: int,
    rows: list[tuple[str, str, str | None]],
) -> None:
    """Assert feature columns section; ignore Rich table padding differences."""
    assert "└── Dataset features" in describe_str
    assert f"└── columns ({n_columns})" in describe_str
    for name, dtype, values in rows:
        if values is None:
            pattern = rf"^\s*{re.escape(name)}\s+{re.escape(dtype)}\s*$"
        else:
            pattern = (
                rf"^\s*{re.escape(name)}\s+{re.escape(dtype)}\s+{re.escape(values)}\s*$"
            )
        assert re.search(pattern, describe_str, re.MULTILINE), (
            f"Missing describe row for {name!r} ({dtype=!r}, {values=!r})"
        )


def test_float_int_casting():
    # this test is only needed for as long as we let JS write data into RecordJson
    # for JS a 3 is a valid float even though any python json parser interprets it as an int
    feature_int = ln.Feature(name="feature_int", dtype=int).save()
    feature_float = ln.Feature(name="feature_float", dtype=float).save()
    test_schema = ln.Schema([feature_int, feature_float], name="test_schema").save()
    sheet = ln.Record(name="TestSheet", is_type=True, schema=test_schema).save()
    record = ln.Record(name="test_record", type=sheet).save()
    record.features.add_values({"feature_int": 5, "feature_float": 3.0})
    record_json = ln.models.RecordJson.get(record=record, feature=feature_float)
    record_json.value = 3
    record_json.save()
    df = sheet.to_dataframe()
    assert df["feature_int"].dtype.name == "int64"
    assert df["feature_float"].dtype.name == "float64"
    # this export call would error if we didn't have type casting
    artifact = sheet.to_artifact()

    related_schemas = list(artifact.schemas.all())
    artifact.schemas.clear()
    artifact.delete(permanent=True)
    record.delete(permanent=True)
    sheet.delete(permanent=True)
    for schema in related_schemas:
        schema.delete(permanent=True)
    # schema.delete(permanent=True), not necessary because already deleted above
    feature_float.delete(permanent=True)
    feature_int.delete(permanent=True)


def test_record_example_compound_treatment(
    populate_sheets_compound_treatment: tuple[ln.Record, ln.Record],  # noqa: F811
):
    treatments_sheet, sample_sheet1 = populate_sheets_compound_treatment

    dictionary = (
        ln.Record.filter(type=treatments_sheet)
        .to_dataframe()[["is_type", "name"]]
        .to_dict(orient="list")
    )
    assert dictionary == {
        "is_type": [
            False,
            False,
        ],
        "name": [
            "treatment2",
            "treatment1",
        ],
    }

    dictionary = (
        ln.Record.filter(type=treatments_sheet)
        .to_dataframe(features=True)[["compound", "concentration", "name"]]
        .to_dict(orient="list")
    )
    assert dictionary == {
        "compound": [
            "drug2",
            "drug1",
        ],
        "concentration": [
            "4nM",
            "2nM",
        ],
        "name": [
            "treatment2",
            "treatment1",
        ],
    }

    partial_df = sample_sheet1.to_dataframe(features=["cell_line", "treatment"])
    assert partial_df.index.name == "name"
    assert partial_df.index.tolist() == ["Sample 1", "Sample 2"]
    dictionary = partial_df[["cell_line", "treatment"]].to_dict(orient="list")
    assert dictionary == {
        "cell_line": [
            "HEK293T",
            "HEK293T",
        ],
        "treatment": [
            "treatment1",
            "treatment2",
        ],
    }

    assert sample_sheet1.input_of_runs.count() == 1
    df = sample_sheet1.to_dataframe()
    assert sample_sheet1.input_of_runs.count() == 2
    assert df.index.name == "name"
    assert df.index.tolist() == ["Sample 1", "Sample 2"]
    assert "name" not in df.columns
    assert not any(col.startswith("__lamindb_record_") for col in df.columns)
    dictionary = df[
        [
            "id",  # a feature
            "uid",  # a feature
            "cell_line",
            "treatment",
            "preparation_date",
        ]
    ].to_dict(orient="list")
    assert dictionary == {
        "id": [1, 2],
        "uid": ["S1", "S2"],
        "cell_line": [
            "HEK293T",
            "HEK293T",
        ],
        "preparation_date": [
            pd.to_datetime("2025-06-01T05:00:00"),
            pd.to_datetime("2025-06-01T06:00:00"),
        ],
        "treatment": [
            "treatment1",
            "treatment2",
        ],
    }

    artifact = sample_sheet1.to_artifact()
    assert sample_sheet1.schema.members.to_list("name") == [
        "name",
        "id",
        "uid",
        "treatment",
        "cell_line",
        "preparation_date",
        "project",
    ]
    assert artifact.run.input_records.count() == 3
    assert artifact.transform.kind == "function"
    assert artifact.transform.key == "__lamindb_record_export__"
    assert artifact.run.status == "completed"
    assert artifact.run.started_at is not None
    assert artifact.run.finished_at is not None
    # looks something like this:
    # name,id,uid,treatment,cell_line,preparation_date,project,__lamindb_record_uid__
    # Sample 1,1,S1,treatment1,HEK293T,2025-06-01 05:00:00,Project 1,iCwgKgZELoLtIoGy
    assert len(artifact.load()) == 2  # two rows in the dataframe
    assert artifact.path.read_text().startswith("""\
name,id,uid,treatment,cell_line,preparation_date,project
Sample 1,1,S1,treatment1,HEK293T,2025-06-01 05:00:00,Project 1""")
    assert artifact.key == f"sheet_exports/{sample_sheet1.name}.csv"
    assert artifact.description.startswith(f"Export of sheet {sample_sheet1.uid}")
    assert artifact._state.adding is False
    assert ln.models.ArtifactRecord.filter(artifact=artifact).count() == 2
    _assert_describe_feature_columns(
        artifact.features.describe(return_str=True),
        7,
        [
            ("cell_line", "bionty.CellLine", "HEK293T"),
            ("id", "int", None),
            ("name", "str", None),
            ("preparation_date", "datetime", None),
            ("project", "Project", "Project 1"),
            ("treatment", "Record[Treatment]", "treatment1, treatment2"),
            ("uid", "str", None),
        ],
    )
    # re-run the export which triggers hash lookup
    artifact2 = sample_sheet1.to_artifact()
    # soft-delete a record in the sheet
    sample_sheet1.records.first().delete()
    assert ln.Record.filter(type=sample_sheet1).count() == 1
    df = sample_sheet1.to_dataframe()
    print(df)
    assert len(df) == 1  # one row in the dataframe

    artifact.delete(permanent=True)
    if artifact2.id != artifact.id:
        artifact2.delete(permanent=True)


def test_nextflow_sheet_with_samples(
    populate_nextflow_sheet_with_samples: ln.Record,  # noqa: F811
):
    """Test the example fixture for nextflow sheet with samples."""
    # This test is to ensure that the fixture works as expected
    # and that the data is correctly populated in the database.
    nextflow_sheet = populate_nextflow_sheet_with_samples

    df = nextflow_sheet.to_dataframe()

    assert df[
        ["expected_cells", "fastq_1", "fastq_2", "sample", "__lamindb_record_name__"]
    ].to_dict(orient="list") == {
        "expected_cells": [
            5000,
            5000,
            5000,
        ],
        "fastq_1": [
            "https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_X_S1_L001_R1_001.fastq.gz",
            "https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_Y_S1_L001_R1_001.fastq.gz",
            "https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_Y_S1_L002_R1_001.fastq.gz",
        ],
        "fastq_2": [
            "https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_X_S1_L001_R2_001.fastq.gz",
            "https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_Y_S1_L001_R2_001.fastq.gz",
            "https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_Y_S1_L002_R2_001.fastq.gz",
        ],
        "__lamindb_record_name__": [
            None,
            None,
            None,
        ],
        "sample": [
            "Sample_X",
            "Sample_Y",
            "Sample_Y",
        ],
    }

    assert nextflow_sheet.schema is not None
    artifact = nextflow_sheet.to_artifact()
    assert artifact.schema is nextflow_sheet.schema
    assert artifact._state.adding is False
    assert set(nextflow_sheet.schema.members.to_list("name")) == {
        "sample",
        "fastq_1",
        "fastq_2",
        "expected_cells",
        "seq_center",
    }
    assert set(artifact.features.slots["columns"].members.to_list("name")) == {
        "sample",
        "fastq_1",
        "fastq_2",
        "expected_cells",
        "seq_center",
    }
    assert artifact.path.read_text().startswith("""\
sample,fastq_1,fastq_2,expected_cells,seq_center,__lamindb_record_uid__,__lamindb_record_name__
Sample_X,https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_X_S1_L001_R1_001.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_X_S1_L001_R2_001.fastq.gz,5000,,""")
    _assert_describe_feature_columns(
        artifact.features.describe(return_str=True),
        5,
        [
            ("expected_cells", "int", None),
            ("fastq_1", "str", None),
            ("fastq_2", "str", None),
            ("sample", "Record[BioSample]", "Sample_X, Sample_Y"),
            ("seq_center", "str", None),
        ],
    )

    related_schemas = list(artifact.schemas.all())
    artifact.schemas.clear()
    for schema in related_schemas:
        schema.delete(permanent=True)
    artifact.delete(permanent=True)


def test_record_soft_deleted_recreate():
    """Test that a soft-deleted record can be recreated with changes."""
    # testing soft delete and recreate with sqlite (postgres is tested in core/test_delete.py)
    # soft delete a record, then recreate it with some changes
    record = ln.Record(name="test_record").save()
    uid = record.uid
    assert record.branch_id == 1
    record.delete()
    assert record.branch_id == -1
    # now recreate the same record with the same uid but a different name
    record = ln.Record(name="test_record 2")
    record.uid = uid
    record.save()
    # now this record is recovered from the trash
    assert record.branch_id == 1
    assert record.name == "test_record 2"
    ln.Record.objects.filter().delete()


def test_annotate_with_user_feature():
    """Test that annotating with a user feature works as expected."""
    user_feature = ln.Feature(name="created_by", dtype=ln.User).save()
    schema = ln.Schema(
        name="test_schema_user_feature",
        features=[user_feature],
        coerce=True,
    ).save()
    sheet = ln.Record(name="A sheet with users", is_type=True, schema=schema).save()
    record = ln.Record(name="first user", type=sheet).save()
    user = ln.User(uid="abcdefgh", handle="test-user").save()
    ln.models.RecordUser(record=record, feature=user_feature, value=user).save()

    df = sheet.to_dataframe()
    assert df.index.name == "__lamindb_record_id__"
    assert df.columns.to_list() == [
        "created_by",
        "__lamindb_record_uid__",
        "__lamindb_record_name__",
    ]
    assert df.iloc[0]["created_by"] == "test-user"

    # clean up
    record.type = None
    record.save()
    record.delete(permanent=True)
    sheet.delete(permanent=True)
    schema.delete(permanent=True)
    user_feature.delete(permanent=True)
    user.delete(permanent=True)


def test_to_artifact_exports_all_records():
    # create sheet with >100 records, the default limit for to_dataframe
    sheet = ln.Record(name="LargeSheet", is_type=True).save()
    for i in range(101):
        ln.Record(name=f"record_{i}", type=sheet).save()
    df = sheet.to_dataframe()
    assert len(df) == 101, f"Expected 101 records, got {len(df)}"
    sheet.records.all().delete(permanent=True)
    sheet.delete(permanent=True)


def test_to_artifact_with_required_non_nullable_data_id_maximal_set_true():
    feature_data_id = ln.Feature(name="data_id", dtype=str, nullable=False).save()
    schema = ln.Schema(
        [feature_data_id],
        name="schema_with_required_data_id",
        maximal_set=True,
    ).save()
    sheet = ln.Record(name="SheetWithDataId", is_type=True, schema=schema).save()
    # Name is intentionally omitted to mirror sheet records in real-world pipelines.
    record = ln.Record(type=sheet).save()
    record.features.add_values({"data_id": "D1"})

    artifact = sheet.to_artifact()
    assert artifact.space_id == sheet.space_id
    assert artifact.run.space_id == sheet.space_id

    df = artifact.load()
    assert "data_id" in df.columns
    assert df["data_id"].to_list() == ["D1"]
    assert "__lamindb_record_name__" in df.columns
    assert df["__lamindb_record_name__"].isna().all()

    # clean up
    export_run = artifact.run
    record.delete(permanent=True)
    sheet.delete(permanent=True)
    artifact.delete(permanent=True)
    export_run.delete(permanent=True)
    schema.delete(permanent=True)
    feature_data_id.delete(permanent=True)


def test_record_export_links_all_upstream_records():
    sheet = ln.Record(name="Run-linked sheet", is_type=True).save()
    record_a = ln.Record(name="record_a", type=sheet).save()
    record_b = ln.Record(name="record_b", type=sheet).save()
    record_c = ln.Record(name="record_c", type=sheet).save()

    transform_a = ln.Transform(key="record_export_upstream_a", kind="function").save()
    transform_b = ln.Transform(key="record_export_upstream_b", kind="function").save()
    run_a = ln.Run(transform=transform_a, status="started").save()
    run_b = ln.Run(transform=transform_b, status="started").save()

    record_a.run = run_a
    record_a.save()
    record_b.run = run_a
    record_b.save()
    record_c.run = run_b
    record_c.save()

    artifact = sheet.to_artifact()
    try:
        linked_record_ids = set(artifact.run.input_records.to_list("id"))
        assert linked_record_ids == {sheet.id, record_a.id, record_b.id, record_c.id}
        linked_run_ids = {
            record.run_id
            for record in artifact.run.input_records.all()
            if record.run_id is not None
        }
        assert linked_run_ids == {run_a.id, run_b.id}
    finally:
        artifact.delete(permanent=True)
        record_a.delete(permanent=True)
        record_b.delete(permanent=True)
        record_c.delete(permanent=True)
        sheet.delete(permanent=True)
        transform_a.delete(permanent=True)
        transform_b.delete(permanent=True)


def test_record_export_links_record_type_when_link_records_false(
    populate_sheets_compound_treatment: tuple[ln.Record, ln.Record],  # noqa: F811
):
    _, sample_sheet = populate_sheets_compound_treatment

    sample_sheet.to_dataframe(link_individual_inputs=False)
    dataframe_export_run = sample_sheet.input_of_runs.order_by("-created_at").first()
    assert dataframe_export_run.input_records.count() == 1
    assert dataframe_export_run.input_records.get().id == sample_sheet.id

    artifact = sample_sheet.to_artifact(link_individual_inputs=False)
    try:
        assert artifact.run.input_records.count() == 1
        assert artifact.run.input_records.get().id == sample_sheet.id
    finally:
        artifact.delete(permanent=True)


def test_record_export_applies_filters():
    sample_sheet = ln.Record(name="FilterSheet", is_type=True).save()
    sample1 = ln.Record(name="sample1", type=sample_sheet).save()
    sample2 = ln.Record(name="sample2", type=sample_sheet).save()

    filtered_df = sample_sheet.to_dataframe(filters={"name": "sample1"})
    assert len(filtered_df) == 1
    assert filtered_df["__lamindb_record_name__"].to_list() == ["sample1"]

    dataframe_export_run = sample_sheet._export_run
    assert dataframe_export_run is not None
    assert dataframe_export_run.input_records.count() == 2
    assert set(dataframe_export_run.input_records.to_list("name")) == {
        "sample1",
        "FilterSheet",
    }

    artifact = sample_sheet.to_artifact(filters={"name": "sample1"})
    try:
        assert len(artifact.load()) == 1
        assert artifact.run is not None
        assert artifact.run.input_records.count() == 2
        assert set(artifact.run.input_records.to_list("name")) == {
            "sample1",
            "FilterSheet",
        }
    finally:
        artifact.delete(permanent=True)
        sample1.delete(permanent=True)
        sample2.delete(permanent=True)
        sample_sheet.delete(permanent=True)


def test_record_export_applies_feature_predicate_filters():
    sample_sheet = ln.Record(name="PredicateFilterSheet", is_type=True).save()
    sample1 = ln.Record(name="sample1", type=sample_sheet).save()
    sample2 = ln.Record(name="sample2", type=sample_sheet).save()
    export_filter_score = ln.Feature(name="export_filter_score", dtype=int).save()
    sample1.features.add_values({"export_filter_score": 10})
    sample2.features.add_values({"export_filter_score": 20})

    filtered_df = sample_sheet.to_dataframe(filters=export_filter_score > 15)
    assert len(filtered_df) == 1
    assert filtered_df["__lamindb_record_name__"].to_list() == ["sample2"]

    dataframe_export_run = sample_sheet._export_run
    assert dataframe_export_run is not None
    assert dataframe_export_run.input_records.count() == 2
    assert set(dataframe_export_run.input_records.to_list("name")) == {
        "sample2",
        "PredicateFilterSheet",
    }

    artifact = sample_sheet.to_artifact(filters=export_filter_score > 15)
    try:
        assert len(artifact.load()) == 1
        assert artifact.run is not None
        assert artifact.run.input_records.count() == 2
        assert set(artifact.run.input_records.to_list("name")) == {
            "sample2",
            "PredicateFilterSheet",
        }
    finally:
        artifact.delete(permanent=True)
        sample1.delete(permanent=True)
        sample2.delete(permanent=True)
        sample_sheet.delete(permanent=True)
        export_filter_score.delete(permanent=True)


def test_record_export_reuses_legacy_transform_uid(
    populate_sheets_compound_treatment: tuple[ln.Record, ln.Record],  # noqa: F811
):
    _, sample_sheet = populate_sheets_compound_treatment
    legacy_transform = ln.Transform(
        key="__lamindb_record_export__",
        kind="function",
    )
    legacy_transform.uid = "LeGaCyUid1230000"
    legacy_transform = legacy_transform.save()

    artifact = sample_sheet.to_artifact()
    try:
        legacy_transform_reloaded = ln.Transform.get(id=legacy_transform.id)
        assert artifact.transform.uid == "v6KpQx9mRt2B0000"
        assert artifact.transform.id == legacy_transform.id
        assert artifact.run is not None
        assert artifact.run.finished_at is not None
        assert artifact.run.status == "completed"
        assert artifact.run.started_at is not None
        assert legacy_transform_reloaded.uid == "v6KpQx9mRt2B0000"
        assert (
            ln.Transform.filter(
                key="__lamindb_record_export__", kind="function"
            ).count()
            == 1
        )
    finally:
        artifact.delete(permanent=True)


def test_record_export_populates_initiated_by_run(
    populate_sheets_compound_treatment: tuple[ln.Record, ln.Record],  # noqa: F811
):
    _, sample_sheet = populate_sheets_compound_treatment
    transform = ln.Transform(key="test_record_export_initiator", kind="function").save()
    ln.track(transform=transform)
    initiating_run = ln.context.run

    artifact = sample_sheet.to_artifact()
    try:
        assert artifact.transform.key == "__lamindb_record_export__"
        assert artifact.run is not None
        assert artifact.run.initiated_by_run is not None
        assert artifact.run.initiated_by_run.id == initiating_run.id
    finally:
        ln.context._run = None
        artifact.delete(permanent=True)


def test_to_artifact_with_duplicate_linked_record_names():
    """Sheet export round-trip resolves Record links by row uid, not ambiguous names."""
    pooled_sample_type = ln.Record(name="PooledSampleDupNameTest", is_type=True).save()
    sample_a = ln.Record(name="poolsample1", type=pooled_sample_type).save()
    sample_b = ln.Record(type=pooled_sample_type).save()
    sample_b.name = "poolsample1"
    sample_b.save()
    assert sample_a.id != sample_b.id

    sample_feature = ln.Feature(
        name="dup_name_pooled_sample", dtype=pooled_sample_type
    ).save()
    fastq_feature = ln.Feature(name="dup_name_fastq_1", dtype=str).save()
    schema = ln.Schema(
        name="duplicate sample names schema",
        features=[sample_feature, fastq_feature],
    ).save()
    sheet_type = ln.Record(name="DuplicateNameSheetType", is_type=True).save()
    sheet = ln.Record(
        name="duplicate name sheet",
        is_type=True,
        type=sheet_type,
        schema=schema,
    ).save()

    row_a = ln.Record(type=sheet).save()
    row_b = ln.Record(type=sheet).save()
    ln.models.RecordRecord(record=row_a, feature=sample_feature, value=sample_a).save()
    ln.models.RecordRecord(record=row_b, feature=sample_feature, value=sample_b).save()
    ln.models.RecordJson(record=row_a, feature=fastq_feature, value="read_a").save()
    ln.models.RecordJson(record=row_b, feature=fastq_feature, value="read_b").save()

    df = sheet.to_dataframe()
    assert df["dup_name_pooled_sample"].tolist() == ["poolsample1", "poolsample1"]
    assert set(df["__lamindb_record_uid__"]) == {row_a.uid, row_b.uid}

    artifact = sheet.to_artifact()
    try:
        assert artifact.schema is sheet.schema
        linked_sample_uids = {
            record.uid
            for record in artifact.records.all()
            if record.type_id == pooled_sample_type.id
        }
        assert linked_sample_uids == {sample_a.uid, sample_b.uid}
    finally:
        export_run = artifact.run
        artifact.delete(permanent=True)
        if export_run is not None:
            export_run.delete(permanent=True)
        row_a.delete(permanent=True)
        row_b.delete(permanent=True)
        sheet.delete(permanent=True)
        sheet_type.delete(permanent=True)
        schema.delete(permanent=True)
        sample_feature.delete(permanent=True)
        fastq_feature.delete(permanent=True)
        sample_a.delete(permanent=True)
        sample_b.delete(permanent=True)
        pooled_sample_type.delete(permanent=True)
