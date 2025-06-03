import lamindb as ln


def test_record_example_compound_treatment():
    structure = ln.Feature(name="structure", dtype="str").save()
    compound_type = ln.Record(name="Compound", is_type=True).save()

    # populate drug1
    drug1 = ln.Record(name="drug1", type=compound_type).save()
    ln.models.RecordJson(record=drug1, feature=structure, value="12345").save()
    # populate drug2
    drug2 = ln.Record(name="drug2", type=compound_type).save()
    ln.models.RecordJson(record=drug2, feature=structure, value="45678").save()

    compound = ln.Feature(name="compound", dtype=compound_type).save()
    concentration = ln.Feature(name="concentration", dtype=compound_type).save()
    sheet1 = ln.Sheet(name="My sheet 1").save()

    treatment_type = ln.Record(name="Treatment", is_type=True).save()
    # populate treatment1
    treatment1 = ln.Record(name="drug1", type=treatment_type, sheet=sheet1).save()
    ln.models.RecordRecord(record=treatment1, feature=compound, value=drug1).save()
    assert drug1 in treatment1.components.all()
    assert treatment1 in drug1.composites.all()
    ln.models.RecordJson(record=treatment1, feature=concentration, value="2nM").save()
    # populate treatment2
    treatment2 = ln.Record(name="drug2", type=treatment_type, sheet=sheet1).save()
    ln.models.RecordRecord(record=treatment2, feature=compound, value=drug2).save()
    ln.models.RecordJson(record=treatment2, feature=concentration, value="4nM").save()
