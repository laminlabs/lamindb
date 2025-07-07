import lamindb as ln
import pandas as pd
from lamindb.examples.fixtures.sheets import (
    populate_nextflow_sheet_with_samples,  # noqa: F401
    populate_sheets_compound_treatment,  # noqa: F401
)


def test_record_example_compound_treatment(
    populate_sheets_compound_treatment: tuple[ln.Record, ln.Record],  # noqa: F811
):
    treatments_sheet, sample_sheet1 = populate_sheets_compound_treatment

    dictionary = (
        ln.Record.filter(type=treatments_sheet)
        .df()[["is_type", "name"]]
        .to_dict(orient="list")
    )
    assert dictionary == {
        "is_type": [
            False,
            False,
        ],
        "name": [
            "treatment1",
            "treatment2",
        ],
    }

    dictionary = (
        ln.Record.filter(type=treatments_sheet)
        .df(features=True)[["compound", "concentration", "name"]]
        .to_dict(orient="list")
    )
    assert dictionary == {
        "compound": [
            "drug1",
            "drug2",
        ],
        "concentration": [
            "2nM",
            "4nM",
        ],
        "name": [
            "treatment1",
            "treatment2",
        ],
    }

    dictionary = (
        ln.Record.filter(type=sample_sheet1)
        .df(features=["cell_line", "treatment"])[["cell_line", "name", "treatment"]]
        .to_dict(orient="list")
    )
    assert dictionary == {
        "cell_line": [
            "HEK293T cell",
            "HEK293T cell",
        ],
        "name": [
            "sample1",
            "sample2",
        ],
        "treatment": [
            "treatment1",
            "treatment2",
        ],
    }

    df = sample_sheet1.to_pandas()
    dictionary = df[
        ["cell_line", "__lamindb_record_name__", "treatment", "preparation_date"]
    ].to_dict(orient="list")
    assert dictionary == {
        "cell_line": [
            "HEK293T cell",
            "HEK293T cell",
        ],
        "__lamindb_record_name__": [
            "sample1",
            "sample2",
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

    # this sheet does not have a schema!
    artifact = sample_sheet1.to_artifact()
    assert sample_sheet1.schema.members.list("name") == [
        "treatment",
        "cell_line",
        "preparation_date",
    ]
    assert artifact.run.input_records.count() == 1
    assert artifact.transform.type == "function"
    # looks something like this:
    # treatment,cell_line,preparation_date,__lamindb_record_uid__,__lamindb_record_name__
    # treatment1,HEK293T cell,2025-06-01 05:00:00,iCwgKgZELoLtIoGy,sample1
    # treatment2,HEK293T cell,2025-06-01 06:00:00,qvU9m7VF6fSdsqJs,sample2
    assert artifact.path.read_text().startswith("""\
treatment,cell_line,preparation_date,__lamindb_record_uid__,__lamindb_record_name__
treatment1,HEK293T cell,2025-06-01 05:00:00""")
    assert artifact.key == f"sheet_exports/{sample_sheet1.name}.csv"
    assert artifact.description.startswith(f"Export of sheet {sample_sheet1.uid}")
    assert artifact._state.adding is False
    assert ln.models.ArtifactRecord.filter(artifact=artifact).count() == 2
    assert (
        artifact.features.describe(return_str=True)
        == """\
Artifact .csv · DataFrame · dataset
└── Dataset features
    └── columns • 3         [Feature]
        cell_line           cat[bionty.CellLine]    HEK293T cell
        treatment           cat[Record[Treatment]]  treatment1, treatment2
        preparation_date    datetime"""
    )
    # re-run the export which triggers hash lookup, which need to escapte re-validation
    sample_sheet1.to_artifact()
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
    assert nextflow_sheet.schema.members.list("name") == [
        "sample",
        "fastq_1",
        "fastq_2",
        "expected_cells",
        "seq_center",
    ]
    assert artifact.features.slots["columns"].members.list("name") == [
        "sample",
        "fastq_1",
        "fastq_2",
        "expected_cells",
    ]
    print(artifact.path.read_text())
    print(artifact.features.describe(return_str=True))
    assert artifact.path.read_text().startswith("""\
sample,fastq_1,fastq_2,expected_cells,__lamindb_record_uid__,__lamindb_record_name__
Sample_X,https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_X_S1_L001_R1_001.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/Sample_X_S1_L001_R2_001.fastq.gz,5000,""")
    assert (
        artifact.features.describe(return_str=True)
        == """\
Artifact .csv · DataFrame · dataset
└── Dataset features
    └── columns • 4         [Feature]
        sample              cat[Record[BioSample]]  Sample_X, Sample_Y
        fastq_1             str
        fastq_2             str
        expected_cells      int"""
    )
    artifact.delete(permanent=True)
