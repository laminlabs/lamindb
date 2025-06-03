import bionty as bt
import lamindb as ln


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
    concentration = ln.Feature(name="concentration", dtype=compound_type).save()
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
    perturbation = ln.Feature(name="perturbation", dtype=compound_type).save()
    cell_line = ln.Feature(name="cell_line", dtype=bt.CellLine).save()
    schema = ln.Schema(
        name="My samples schema 2025-05", features=[perturbation, cell_line]
    ).save()
    # a sheet for samples
    sheet1 = ln.Sheet(name="My samples 2025-05", schema=schema).save()
    # values for cell lines
    hek293t = bt.CellLine.from_source(name="HEK293T").save()

    # populate sample1
    sample1 = ln.Record(name="sample1", type=sample_type, sheet=sheet1).save()
    ln.models.RecordRecord(
        record=sample1, feature=perturbation, value=treatment1
    ).save()
    bt.models.RecordCellLine(
        record=sample1, feature=cell_line, cellline=hek293t
    ).save()  # parallel to ArtifactCellLine
    # populate sample2
    sample2 = ln.Record(name="sample2", type=sample_type, sheet=sheet1).save()
    ln.models.RecordRecord(
        record=sample2, feature=perturbation, value=treatment2
    ).save()
    bt.models.RecordCellLine(
        record=sample2, feature=cell_line, cellline=hek293t
    ).save()  # parallel to ArtifactCellLine

    # another sheet for samples
    sample_note = ln.Feature(name="sample_note", dtype="str").save()
    schema2 = ln.Schema(
        name="My samples schema 2025-06",
        features=[perturbation, cell_line, sample_note],
    ).save()
    # the sheet
    sheet2 = ln.Sheet(name="My samples 2025-06", schema=schema2).save()
    # populate sample3
    sample3 = ln.Record(type=sample_type, sheet=sheet2).save()  # no name
    ln.models.RecordRecord(
        record=sample3, feature=perturbation, value=treatment1
    ).save()
    bt.models.RecordCellLine(
        record=sample3, feature=cell_line, cellline=hek293t
    ).save()  # parallel to ArtifactCellLine
    # populate sample4
    sample4 = ln.Record(type=sample_type, sheet=sheet2).save()
    ln.models.RecordRecord(
        record=sample4, feature=perturbation, value=treatment2
    ).save()
    bt.models.RecordCellLine(
        record=sample4, feature=cell_line, cellline=hek293t
    ).save()  # parallel to ArtifactCellLine
