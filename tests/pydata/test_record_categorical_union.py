import lamindb as ln
from django.db import transaction


def test_record_to_dataframe_categorical_union():
    with transaction.atomic():
        label = ln.ULabel(name="union_export_label").save()
        project = ln.Project(name="union_export_project").save()
        feature = ln.Feature(
            name="union_export_category", dtype="cat[ULabel.name|Project.name]"
        ).save()
        schema = ln.Schema([feature]).save()
        sheet = ln.Record(name="UnionExportSheet", is_type=True, schema=schema).save()
        expected = {"label-row": label.name, "project-row": project.name}
        records = [
            ln.Record(name=name, type=sheet, features={feature.name: value}).save()
            for name, value in expected.items()
        ]

        assert {
            record.name: record.features.get_values()[feature.name]
            for record in records
        } == expected

        exported = sheet.to_dataframe(is_run_input=False).set_index(
            "__lamindb_record_name__"
        )
        assert exported[feature.name].to_dict() == expected

        transaction.set_rollback(True)
