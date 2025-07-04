import lamindb as ln
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

    dictionary = (
        ln.Record.filter(type=sample_sheet1)
        .df(features="queryset")[["cell_line", "name", "treatment", "preparation_date"]]
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
        "preparation_date": [
            {
                "2025-06-01T05:00:00Z",
            },
            {
                "2025-06-01T06:00:00Z",
            },
        ],
        "treatment": [
            "treatment1",
            "treatment2",
        ],
    }


def test_nextflow_sheet_with_samples(
    populate_nextflow_sheet_with_samples: ln.Record,  # noqa: F811
):
    """Test the example fixture for nextflow sheet with samples."""
    # This test is to ensure that the fixture works as expected
    # and that the data is correctly populated in the database.
    nextflow_sheet = populate_nextflow_sheet_with_samples

    assert ln.Record.filter(type=nextflow_sheet).df(features="queryset")[
        ["expected_cells", "fastq_1", "fastq_2", "sample", "name"]
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
        "name": [
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
