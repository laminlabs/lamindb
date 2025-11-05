import lamindb as ln
import pandas as pd
from lamindb.examples.fixtures.sheets import (
    populate_nextflow_sheet_with_samples,  # noqa: F401
    populate_sheets_compound_treatment,  # noqa: F401
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
    df = sheet.type_to_dataframe()
    assert df["feature_int"].dtype.name == "int64"
    assert df["feature_float"].dtype.name == "float64"
    # this export call would error if we didn't have type casting
    artifact = sheet.to_artifact()

    related_schemas = list(artifact.feature_sets.all())
    artifact.feature_sets.clear()
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

    dictionary = (
        ln.Record.filter(type=sample_sheet1)
        .to_dataframe(features=["cell_line", "treatment"])[
            ["cell_line", "__lamindb_record_name__", "treatment"]
        ]
        .to_dict(orient="list")
    )
    assert dictionary == {
        "cell_line": [
            "HEK293T cell",
            "HEK293T cell",
        ],
        "__lamindb_record_name__": [
            "sample2",
            "sample1",
        ],
        "treatment": [
            "treatment2",
            "treatment1",
        ],
    }

    df = sample_sheet1.type_to_dataframe()
    assert df.index.name == "__lamindb_record_id__"
    dictionary = df[
        [
            "id",  # a feature
            "uid",  # a feature
            "name",  # a feature
            "cell_line",
            "treatment",
            "preparation_date",
            "__lamindb_record_name__",
        ]
    ].to_dict(orient="list")
    assert dictionary == {
        "id": [1, 2],
        "uid": ["S1", "S2"],
        "name": ["Sample 1", "Sample 2"],
        "cell_line": [
            "HEK293T cell",
            "HEK293T cell",
        ],
        "preparation_date": [
            pd.to_datetime("2025-06-01T05:00:00"),
            pd.to_datetime("2025-06-01T06:00:00"),
        ],
        "treatment": [
            "treatment1",
            "treatment2",
        ],
        "__lamindb_record_name__": [
            "sample1",
            "sample2",
        ],
    }

    artifact = sample_sheet1.to_artifact()
    assert sample_sheet1.schema.members.to_list("name") == [
        "id",
        "uid",
        "name",
        "treatment",
        "cell_line",
        "preparation_date",
        "project",
    ]
    assert artifact.run.input_records.count() == 1
    assert artifact.transform.type == "function"
    # looks something like this:
    # id,uid,name,treatment,cell_line,preparation_date,__lamindb_record_uid__,__lamindb_record_name__
    # 1,S1,Sample 1,treatment1,HEK293T cell,2025-06-01 05:00:00,iCwgKgZELoLtIoGy,sample1
    # 2,S2,Sample 2,treatment2,HEK293T cell,2025-06-01 06:00:00,qvU9m7VF6fSdsqJs,sample2
    assert len(artifact.load()) == 2  # two rows in the dataframe
    assert artifact.path.read_text().startswith("""\
id,uid,name,treatment,cell_line,preparation_date,project,__lamindb_record_uid__,__lamindb_record_name__
1,S1,Sample 1,treatment1,HEK293T cell,2025-06-01 05:00:00,Project 1""")
    assert artifact.key == f"sheet_exports/{sample_sheet1.name}.csv"
    assert artifact.description.startswith(f"Export of sheet {sample_sheet1.uid}")
    assert artifact._state.adding is False
    assert ln.models.ArtifactRecord.filter(artifact=artifact).count() == 2
    assert artifact.features.describe(return_str=True).endswith("""\
└── Dataset features
    └── columns (7)
        cell_line           bionty.CellLine         HEK293T cell
        project             Project                 Project 1
        treatment           Record[Treatment]       treatment1, treatment2
        id                  int
        uid                 str
        name                str
        preparation_date    datetime""")
    # re-run the export which triggers hash lookup
    sample_sheet1.to_artifact()
    # soft-delete a record in the sheet
    sample_sheet1.records.first().delete()
    assert ln.Record.filter(type=sample_sheet1).count() == 1
    df = sample_sheet1.type_to_dataframe()
    print(df)
    assert len(df) == 1  # one row in the dataframe

    artifact.delete(permanent=True)


def test_nextflow_sheet_with_samples(
    populate_nextflow_sheet_with_samples: ln.Record,  # noqa: F811
):
    """Test the example fixture for nextflow sheet with samples."""
    # This test is to ensure that the fixture works as expected
    # and that the data is correctly populated in the database.
    nextflow_sheet = populate_nextflow_sheet_with_samples

    df = nextflow_sheet.to_pandas()

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
    assert nextflow_sheet.schema.members.to_list("name") == [
        "sample",
        "fastq_1",
        "fastq_2",
        "expected_cells",
        "seq_center",
    ]
    assert artifact.features.slots["columns"].members.to_list("name") == [
        "sample",
        "fastq_1",
        "fastq_2",
        "expected_cells",
    ]
    assert artifact.path.read_text().startswith("""\
sample,fastq_1,fastq_2,expected_cells,__lamindb_record_uid__,__lamindb_record_name__
Sample_X,https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_X_S1_L001_R1_001.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_X_S1_L001_R2_001.fastq.gz,5000,""")
    assert artifact.features.describe(return_str=True).endswith("""\
└── Dataset features
    └── columns (4)
        sample              Record[BioSample]       Sample_X, Sample_Y
        fastq_1             str
        fastq_2             str
        expected_cells      int""")

    related_schemas = list(artifact.feature_sets.all())
    artifact.feature_sets.clear()
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
        coerce_dtype=True,
    ).save()
    sheet = ln.Record(name="A sheet with users", is_type=True, schema=schema).save()
    record = ln.Record(name="first user", type=sheet).save()
    user = ln.User(uid="abcdefgh", handle="test-user").save()
    ln.models.RecordUser(record=record, feature=user_feature, value=user).save()

    df = sheet.type_to_dataframe()
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
