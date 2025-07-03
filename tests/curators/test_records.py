import bionty as bt
import lamindb as ln
import pandas as pd


def test_record_example_compound_treatment(populate_sheets_compound_treatment):
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


def test_record_nextflow_samples():
    # Biosample schema and type
    samples_schema = ln.Schema(
        name="Biosample test schema",
        features=[
            ln.Feature(name="species", dtype="cat[bionty.Organism]").save(),
            ln.Feature(name="cell_type", dtype="cat[bionty.CellType]").save(),
            ln.Feature(name="tissue", dtype="cat[bionty.Tissue]").save(),
        ],
    ).save()

    biosample_type = ln.Record(name="BioSample", is_type=True).save()

    # Biosamples sheet
    samples_sheet = ln.Record(
        name="My samples 2025-04", schema=samples_schema, type=biosample_type
    ).save()
    sample_x = ln.Record(name="Sample_X", type=samples_sheet).save()
    sample_y = ln.Record(name="Sample_Y", type=samples_sheet).save()

    organism_human = bt.Organism.from_source(name="human").save()
    celltype_tcell = bt.CellType.from_source(name="T cell").save()
    tissue_blood = bt.Tissue.from_source(name="blood").save()

    features = ln.Feature.lookup()
    for sample in [sample_x, sample_y]:
        bt.models.RecordOrganism(
            record=sample, feature=features.species, organism=organism_human
        ).save()
        bt.models.RecordCellType(
            record=sample, feature=features.cell_type, celltype=celltype_tcell
        ).save()
        bt.models.RecordTissue(
            record=sample, feature=features.tissue, tissue=tissue_blood
        ).save()

    # Nextflow samplesheet schema
    nextflow_schema = ln.Schema(
        name="RNA-seq standard",
        features=[
            ln.Feature(name="sample", dtype="cat[Record[BioSample]]").save(),
            ln.Feature(name="fastq_1", dtype=str).save(),
            ln.Feature(name="fastq_2", dtype=str).save(),
            ln.Feature(name="expected_cells", dtype=int).save(),
            ln.Feature(name="seq_center", dtype=str).save(),
        ],
    ).save()

    nextflowsample_type = ln.Record(name="NextflowSample", is_type=True).save()
    nextflow_sheet = ln.Record(
        schema=nextflow_schema,
        name="RNA-seq nextflow samplesheet 001",
        type=nextflowsample_type,
        is_type=True,
    ).save()

    df = pd.DataFrame(
        {
            "sample": ["Sample_X", "Sample_Y", "Sample_Y"],
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
            "expected_cells": [5000, 5000, 5000],
        }
    )

    features = ln.Feature.lookup()
    for _, row in df.iterrows():
        sample = ln.Record(type=nextflow_sheet).save()
        ln.models.RecordRecord(
            record=sample,
            feature=features.sample,
            value=ln.Record.get(name=row["sample"]),
        ).save()
        ln.models.RecordJson(
            record=sample, feature=features.fastq_1, value=row["fastq_1"]
        ).save()
        ln.models.RecordJson(
            record=sample, feature=features.fastq_2, value=row["fastq_2"]
        ).save()
        ln.models.RecordJson(
            record=sample, feature=features.expected_cells, value=row["expected_cells"]
        ).save()
