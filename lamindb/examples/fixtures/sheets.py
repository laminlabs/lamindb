import bionty as bt
import pandas as pd
import pytest

import lamindb as ln


@pytest.fixture(scope="module")
def populate_sheets_compound_treatment():
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
    treatments_sheet = ln.Record(
        name="My treatments 2025-05", type=treatment_type, is_type=True
    ).save()  # sheet without validating schema

    # populate treatment1
    treatment1 = ln.Record(name="treatment1", type=treatments_sheet).save()
    ln.models.RecordRecord(record=treatment1, feature=compound, value=drug1).save()
    assert drug1 in treatment1.components.all()  # noqa: S101
    assert treatment1 in drug1.composites.all()  # noqa: S101
    ln.models.RecordJson(record=treatment1, feature=concentration, value="2nM").save()
    # populate treatment2
    treatment2 = ln.Record(name="treatment2", type=treatments_sheet).save()
    ln.models.RecordRecord(record=treatment2, feature=compound, value=drug2).save()
    ln.models.RecordJson(record=treatment2, feature=concentration, value="4nM").save()

    # Samples ---------------------------

    # features named id, uid or name conflict with django field names, we test them here
    id_feature = ln.Feature(name="id", dtype=int).save()
    uid_feature = ln.Feature(name="uid", dtype=str).save()
    name_feature = ln.Feature(name="name", dtype=str).save()

    project = ln.Feature(name="project", dtype=ln.Project).save()
    project1 = ln.Project(name="Project 1").save()
    sample_type = ln.Record(name="BioSample", is_type=True).save()
    treatment = ln.Feature(name="treatment", dtype=treatment_type).save()
    cell_line = ln.Feature(name="cell_line", dtype=bt.CellLine).save()
    preparation_date = ln.Feature(name="preparation_date", dtype="datetime").save()
    cell_line.dtype = "cat[bionty.CellLine]"  # might have previously been set to "cat"
    cell_line.save()
    sample_schema1 = ln.Schema(
        name="My samples schema 2025-06",
        features=[
            id_feature,
            uid_feature,
            name_feature,
            treatment,
            cell_line,
            preparation_date,
            project,
        ],
    ).save()
    sample_sheet1 = ln.Record(
        name="My samples 2025-06", schema=sample_schema1, type=sample_type
    ).save()
    # values for cell lines
    hek293t = bt.CellLine.from_source("HEK293T").save()

    # populate sample1
    sample1 = ln.Record(name="sample1", type=sample_sheet1).save()
    ln.models.RecordJson(record=sample1, feature=id_feature, value=1).save()
    ln.models.RecordJson(record=sample1, feature=uid_feature, value="S1").save()
    ln.models.RecordJson(record=sample1, feature=name_feature, value="Sample 1").save()
    ln.models.RecordRecord(record=sample1, feature=treatment, value=treatment1).save()
    bt.models.RecordCellLine(record=sample1, feature=cell_line, value=hek293t).save()
    ln.models.RecordJson(
        record=sample1, feature=preparation_date, value="2025-06-01T05:00:00"
    ).save()
    ln.models.RecordProject(record=sample1, feature=project, value=project1).save()
    # populate sample2
    sample2 = ln.Record(name="sample2", type=sample_sheet1).save()
    ln.models.RecordJson(record=sample2, feature=id_feature, value=2).save()
    ln.models.RecordJson(record=sample2, feature=uid_feature, value="S2").save()
    ln.models.RecordJson(record=sample2, feature=name_feature, value="Sample 2").save()
    ln.models.RecordRecord(record=sample2, feature=treatment, value=treatment2).save()
    bt.models.RecordCellLine(record=sample2, feature=cell_line, value=hek293t).save()
    ln.models.RecordJson(
        record=sample2, feature=preparation_date, value="2025-06-01T06:00:00"
    ).save()
    ln.models.RecordProject(record=sample2, feature=project, value=project1).save()

    # another sheet for samples
    sample_note = ln.Feature(name="sample_note", dtype="str").save()
    sample_schema2 = ln.Schema(
        name="My samples schema 2025-07",
        features=[treatment, cell_line, sample_note, project],
    ).save()
    # the sheet
    sample_sheet2 = ln.Record(
        name="My samples 2025-07", schema=sample_schema2, type=sample_type
    ).save()
    # populate sample3
    sample3 = ln.Record(type=sample_sheet2).save()  # no name
    ln.models.RecordRecord(record=sample3, feature=treatment, value=treatment1).save()
    bt.models.RecordCellLine(record=sample3, feature=cell_line, value=hek293t).save()
    ln.models.RecordJson(
        record=sample3, feature=preparation_date, value="2025-06-02T05:00:00Z"
    ).save()
    ln.models.RecordProject(record=sample3, feature=project, value=project1).save()
    # populate sample4
    sample4 = ln.Record(type=sample_sheet2).save()
    ln.models.RecordRecord(record=sample4, feature=treatment, value=treatment2).save()
    bt.models.RecordCellLine(record=sample4, feature=cell_line, value=hek293t).save()
    ln.models.RecordJson(
        record=sample4, feature=preparation_date, value="2025-06-02T06:00:00Z"
    ).save()
    ln.models.RecordProject(record=sample4, feature=project, value=project1).save()

    yield treatments_sheet, sample_sheet1

    sample4.delete(permanent=True)
    sample3.delete(permanent=True)
    sample_sheet2.delete(permanent=True)
    sample_schema2.delete(permanent=True)
    sample_note.delete(permanent=True)
    sample2.delete(permanent=True)
    sample1.delete(permanent=True)
    # hek293t.delete(permanent=True)  # not for now
    sample_sheet1.delete(permanent=True)
    sample_schema1.delete(permanent=True)
    preparation_date.delete(permanent=True)
    cell_line.delete(permanent=True)
    # sample_type.delete(permanent=True)   # not for now
    treatment2.delete(permanent=True)
    treatment1.delete(permanent=True)
    treatments_sheet.delete(permanent=True)
    treatment_type.delete(permanent=True)
    concentration.delete(permanent=True)
    drug2.delete(permanent=True)
    drug1.delete(permanent=True)
    structure.delete(permanent=True)
    compound.delete(permanent=True)
    compound_type.delete(permanent=True)


@pytest.fixture(scope="module")
def populate_nextflow_sheet_with_samples():
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
            record=sample, feature=features.species, value=organism_human
        ).save()
        bt.models.RecordCellType(
            record=sample, feature=features.cell_type, value=celltype_tcell
        ).save()
        bt.models.RecordTissue(
            record=sample, feature=features.tissue, value=tissue_blood
        ).save()

    # Nextflow samplesheet schema
    nextflow_schema = ln.Schema(
        name="RNA-seq standard",
        features=[
            ln.Feature(name="sample", dtype="cat[Record[BioSample]]").save(),
            ln.Feature(name="fastq_1", dtype=str).save(),
            ln.Feature(name="fastq_2", dtype=str).save(),
            ln.Feature(name="expected_cells", dtype=int).save(),
            ln.Feature(name="seq_center", dtype=str).save().with_config(optional=True),
        ],
        ordered_set=True,
    ).save()

    nextflowsample_type = ln.Record(name="NextflowSample", is_type=True).save()
    nextflow_sheet = ln.Record(
        schema=nextflow_schema,
        name="RNA-seq nextflow samplesheet 001",
        type=nextflowsample_type,
        is_type=True,
    ).save()

    sample_data = {
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
    df = pd.DataFrame(sample_data)

    features = ln.Feature.lookup()
    nextflow_samples = []
    for _, row in df.iterrows():
        sample = ln.Record(type=nextflow_sheet).save()
        nextflow_samples.append(sample)
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

    yield nextflow_sheet

    # Delete in reverse order of creation
    # Delete nextflow samples
    for sample in reversed(nextflow_samples):
        sample.delete(permanent=True)

    # Delete nextflow sheet and schema
    nextflow_sheet.delete(permanent=True)
    nextflowsample_type.delete(permanent=True)
    nextflow_schema.delete(permanent=True)

    # Delete samples sheet and schema
    samples_sheet.records.all().delete(permanent=True)
    samples_sheet.delete(permanent=True)
    # biosample_type.delete(permanent=True)  # not for now (shared with first fixture)
    samples_schema.delete(permanent=True)

    print(ln.Schema.to_dataframe())

    # Delete nextflow schema features
    features = ln.Feature.lookup()
    features.seq_center.delete(permanent=True)
    features.expected_cells.delete(permanent=True)
    features.fastq_2.delete(permanent=True)
    features.fastq_1.delete(permanent=True)
    features.sample.delete(permanent=True)

    # Delete biosamples
    sample_y.delete(permanent=True)
    sample_x.delete(permanent=True)

    # Delete biosample schema features
    features.tissue.delete(permanent=True)
    features.cell_type.delete(permanent=True)
    features.species.delete(permanent=True)

    # Note: organism_human, celltype_tcell, tissue_blood are from bionty
    # and might be shared, so not deleting them (similar to hek293t in first fixture)
