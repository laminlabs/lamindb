import bionty as bt
import lamindb as ln
import pandas as pd


def test_record_example_compound_treatment():
    # Compounds ---------------------------

    compound_type = ln.Record(name="Compound", is_type=True).save()

    # features for compounds
    structure = ln.Feature(name="structure", dtype="str").save()

    # drug1
    drug1 = ln.Record(name="drug1", type=compound_type).save()
    ln.models.RecordJson(record=drug1, feature=structure, value="12345").save()
    # drug2
    drug2 = ln.Record(name="drug2", type=compound_type).save()
    ln.models.RecordJson(record=drug2, feature=structure, value="45678").save()

    # Treatments ---------------------------

    treatment_type = ln.Record(name="Treatment", is_type=True).save()

    # features for treatments
    compound = ln.Feature(name="compound", dtype=compound_type).save()
    concentration = ln.Feature(name="concentration", dtype="num").save()
    # a sheet for treatments
    my_treatments = ln.Sheet(
        name="My treatments 2025-05"
    ).save()  # sheet without schema

    # populate treatment1
    treatment1 = ln.Record(
        name="drug1", type=treatment_type, sheet=my_treatments
    ).save()
    ln.models.RecordRecord(record=treatment1, feature=compound, value=drug1).save()
    assert drug1 in treatment1.components.all()
    assert treatment1 in drug1.composites.all()
    ln.models.RecordJson(record=treatment1, feature=concentration, value="2nM").save()
    # populate treatment2
    treatment2 = ln.Record(
        name="drug2", type=treatment_type, sheet=my_treatments
    ).save()
    ln.models.RecordRecord(record=treatment2, feature=compound, value=drug2).save()
    ln.models.RecordJson(record=treatment2, feature=concentration, value="4nM").save()

    # Samples ---------------------------

    sample_type = ln.Record(name="BioSample", is_type=True).save()

    # features for samples
    treatment = ln.Feature(name="treatment", dtype=compound_type).save()
    cell_line = ln.Feature(name="cell_line", dtype=bt.CellLine).save()
    schema = ln.Schema(
        name="My samples schema 2025-05", features=[treatment, cell_line]
    ).save()
    # a sheet for samples
    sheet1 = ln.Sheet(name="My samples 2025-05", schema=schema).save()
    # values for cell lines
    hek293t = bt.CellLine.from_source(name="HEK293T").save()

    # populate sample1
    sample1 = ln.Record(name="sample1", type=sample_type, sheet=sheet1).save()
    ln.models.RecordRecord(record=sample1, feature=treatment, value=treatment1).save()
    bt.models.RecordCellLine(
        record=sample1, feature=cell_line, cellline=hek293t
    ).save()  # parallel to ArtifactCellLine
    # populate sample2
    sample2 = ln.Record(name="sample2", type=sample_type, sheet=sheet1).save()
    ln.models.RecordRecord(record=sample2, feature=treatment, value=treatment2).save()
    bt.models.RecordCellLine(
        record=sample2, feature=cell_line, cellline=hek293t
    ).save()  # parallel to ArtifactCellLine

    # another sheet for samples
    sample_note = ln.Feature(name="sample_note", dtype="str").save()
    schema2 = ln.Schema(
        name="My samples schema 2025-06",
        features=[treatment, cell_line, sample_note],
    ).save()
    # the sheet
    sheet2 = ln.Sheet(name="My samples 2025-06", schema=schema2).save()
    # populate sample3
    sample3 = ln.Record(type=sample_type, sheet=sheet2).save()  # no name
    ln.models.RecordRecord(record=sample3, feature=treatment, value=treatment1).save()
    bt.models.RecordCellLine(
        record=sample3, feature=cell_line, cellline=hek293t
    ).save()  # parallel to ArtifactCellLine
    # populate sample4
    sample4 = ln.Record(type=sample_type, sheet=sheet2).save()
    ln.models.RecordRecord(record=sample4, feature=treatment, value=treatment2).save()
    bt.models.RecordCellLine(
        record=sample4, feature=cell_line, cellline=hek293t
    ).save()  # parallel to ArtifactCellLine


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
    samples_sheet = ln.Sheet(name="My samples 2025-04", schema=samples_schema).save()
    sample_x = ln.Record(
        name="Sample_X", type=biosample_type, sheet=samples_sheet
    ).save()
    sample_y = ln.Record(
        name="Sample_Y", type=biosample_type, sheet=samples_sheet
    ).save()

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
    nextflow_sheet = ln.Sheet(
        schema=nextflow_schema, name="RNA-seq nextflow samplesheet 001"
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
        sample = ln.Record(sheet=nextflow_sheet, type=nextflowsample_type).save()
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


if __name__ == "__main__":
    test_record_example_compound_treatment()
