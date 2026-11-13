import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from django.db import transaction


@pytest.fixture(autouse=True)
def isolate_records():
    with transaction.atomic():
        yield
        transaction.set_rollback(True)


@pytest.mark.parametrize(
    "dtype", ["cat[ULabel.name|Project.name]", "cat[Project|ULabel]"]
)
@pytest.mark.parametrize("second_only", [False, True])
def test_record_to_dataframe_categorical_union(dtype, second_only):
    label = ln.ULabel(name="union_export_label").save()
    project = ln.Project(name="union_export_project").save()
    feature = ln.Feature(name="union_export_category", dtype=dtype).save()
    schema = ln.Schema([feature]).save()
    sheet = ln.Record(name="UnionExportSheet", is_type=True, schema=schema).save()
    expected = {"label-row": label.name, "project-row": project.name}
    if second_only:
        expected = {"project-row": project.name}
    records = [
        ln.Record(name=name, type=sheet, features={feature.name: value}).save()
        for name, value in expected.items()
    ]

    assert {
        record.name: record.features.get_values()[feature.name] for record in records
    } == expected

    exported = sheet.to_dataframe(is_run_input=False).set_index(
        "__lamindb_record_name__"
    )
    assert exported[feature.name].to_dict() == expected

    for options, name_column in (
        ({"include": "features"}, "__lamindb_record_name__"),
        ({"features": [feature.name]}, "name"),
    ):
        frame = ln.Record.filter(type=sheet).to_dataframe(**options, is_run_input=False)
        assert frame.set_index(name_column)[feature.name].to_dict() == expected


def test_record_union_export_registry_fields():
    label = ln.ULabel(name="union_field_label").save()
    project = ln.Project(name="union_field_project").save()
    feature = ln.Feature(
        name="union_field", dtype="cat[ULabel.name|Project.uid]"
    ).save()
    sheet = ln.Record(
        name="UnionFieldSheet", is_type=True, schema=ln.Schema([feature]).save()
    ).save()
    expected = {"label": label.name, "project": project.uid}
    for name, value in expected.items():
        ln.Record(name=name, type=sheet, features={feature.name: value}).save()

    frame = sheet.to_dataframe(is_run_input=False)
    assert (
        frame.set_index("__lamindb_record_name__")[feature.name].to_dict() == expected
    )


@pytest.mark.parametrize("list_dtype", [False, True])
def test_record_union_export_merges_three_registries(list_dtype):
    label = ln.ULabel(name="union_merge_label").save()
    project = ln.Project(name="union_merge_project").save()
    other = ln.Record(name="union_merge_record").save()
    dtype = "cat[ULabel.name|Project.name|Record.name]"
    feature = ln.Feature(
        name="union_merge", dtype=f"list[{dtype}]" if list_dtype else dtype
    ).save()
    extra = ln.Feature(name="union_merge_extra", dtype=ln.Project).save()
    schema = ln.Schema([feature, extra]).save()
    sheet = ln.Record(name="UnionMergeSheet", is_type=True, schema=schema).save()
    expected = [label.name, project.name, other.name]
    ln.Record(
        name="all-registries",
        type=sheet,
        features={feature.name: expected, extra.name: project.name},
    ).save()

    frame = sheet.to_dataframe(is_run_input=False)
    actual = frame[feature.name].iloc[0]
    assert isinstance(actual, list if list_dtype else set)
    assert len(actual) == len(expected)
    assert set(actual) == set(expected)
    assert frame[extra.name].iloc[0] == project.name


def test_record_union_export_empty_and_deleted_values():
    label = ln.ULabel(name="union_missing_label").save()
    project = ln.Project(name="union_deleted_project").save()
    feature = ln.Feature(name="union_missing", dtype="cat[ULabel|Project]").save()
    sheet = ln.Record(
        name="UnionMissingSheet", is_type=True, schema=ln.Schema([feature]).save()
    ).save()
    empty = sheet.to_dataframe(is_run_input=False)
    assert empty.empty
    assert isinstance(empty[feature.name].dtype, pd.CategoricalDtype)

    ln.Record(name="missing", type=sheet, features={feature.name: None}).save()
    for name, value in [("label", label.name), ("deleted", project.name)]:
        ln.Record(name=name, type=sheet, features={feature.name: value}).save()
    project.branch_id = -1
    project.save()

    frame = sheet.to_dataframe(is_run_input=False).set_index("__lamindb_record_name__")
    assert frame.loc["label", feature.name] == label.name
    assert pd.isna(frame.loc["missing", feature.name])
    assert pd.isna(frame.loc["deleted", feature.name])
    assert isinstance(frame[feature.name].dtype, pd.CategoricalDtype)


def test_record_union_export_from_dataframe():
    bt.Tissue(name="union tissue", ontology_id="UBERON:0000001").save()
    bt.CellType(name="union cell type", ontology_id="CL:0000904").save()
    feature = ln.Feature(
        name="union_biosample_id",
        dtype="cat[bionty.Tissue.ontology_id|bionty.CellType.ontology_id]",
    ).save()
    sheet = ln.Record(
        name="UnionStudy", is_type=True, schema=ln.Schema([feature]).save()
    ).save()
    expected = {"tissue-study": "UBERON:0000001", "cell-study": "CL:0000904"}
    frame = pd.DataFrame(
        {"study_id": list(expected), feature.name: pd.Categorical(expected.values())}
    )
    ln.Record.from_dataframe(frame, type=sheet, name_field="study_id").save()
    assert {
        record.name: record.features.get_values()[feature.name]
        for record in sheet.records.all()
    } == expected

    exported = sheet.to_dataframe(is_run_input=False).set_index(
        "__lamindb_record_name__"
    )
    assert exported[feature.name].to_dict() == expected


def test_artifact_union_export(tmp_path):
    label = ln.ULabel(name="union_host_label").save()
    project = ln.Project(name="union_host_project").save()
    feature = ln.Feature(name="union_host_value", dtype="cat[ULabel|Project]").save()
    expected = {}
    for value in (label.name, project.name):
        path = tmp_path / f"{value}.txt"
        path.write_text(value)
        artifact = ln.Artifact(path, key=path.name, run=False).save()
        artifact.features.add_values({feature.name: value})
        expected[artifact.uid] = value

    frame = ln.Artifact.filter(uid__in=expected).to_dataframe(features=[feature.name])
    assert frame.set_index("uid")[feature.name].to_dict() == expected


def test_run_union_export(tmp_path):
    record = ln.Record(name="union_run_record").save()
    path = tmp_path / "union_run_artifact.txt"
    path.write_text("union run input")
    artifact = ln.Artifact(path, key=path.name, run=False).save()
    feature = ln.Feature(
        name="union_run_value", dtype="cat[Record.name|Artifact.key]"
    ).save()
    transform = ln.Transform(key="union-export-test").save()
    expected = {}
    for value in (record.name, artifact.key):
        run = ln.Run(transform).save()
        run.features.add_values({feature.name: value})
        expected[run.uid] = value

    frame = ln.Run.filter(uid__in=expected).to_dataframe(features=[feature.name])
    assert frame.set_index("uid")[feature.name].to_dict() == expected


def test_record_union_export_default_field_without_branch():
    label = ln.ULabel(name="union_user_label").save()
    user = ln.User.get(handle=ln.setup.settings.user.handle)
    feature = ln.Feature(name="union_user", dtype="cat[ULabel|User]").save()
    sheet = ln.Record(
        name="UnionUserSheet", is_type=True, schema=ln.Schema([feature]).save()
    ).save()
    expected = {"label": label.name, "user": user.handle}
    for name, value in expected.items():
        ln.Record(name=name, type=sheet, features={feature.name: value}).save()

    frame = sheet.to_dataframe(is_run_input=False)
    assert (
        frame.set_index("__lamindb_record_name__")[feature.name].to_dict() == expected
    )


def test_record_union_export_keeps_selected_feature_id():
    label = ln.ULabel(name="union_duplicate_label").save()
    project = ln.Project(name="union_duplicate_project").save()
    feature = ln.Feature(name="union_duplicate", dtype="cat[ULabel|Project]").save()
    feature_type = ln.Feature(name="UnionDuplicateType", is_type=True).save()
    duplicate = ln.Feature(
        name=feature.name, dtype="cat[ULabel|Project]", type=feature_type
    ).save()
    record = ln.Record(name="union-duplicate-row").save()
    ln.models.RecordULabel(record=record, feature=feature, value=label).save()
    ln.models.RecordProject(record=record, feature=duplicate, value=project).save()

    frame = ln.Record.filter(id=record.id).to_dataframe(
        include="features", is_run_input=False
    )
    assert frame[feature.name].iloc[0] == label.name


@pytest.mark.parametrize("include_feature_id", [False, True])
@pytest.mark.parametrize("union_dtype", [False, True])
def test_record_feature_export_keeps_incoming_links_separate(
    include_feature_id, union_dtype
):
    feature = ln.Feature(
        name="incoming_export_relation",
        dtype="cat[Record.name|ULabel.name]" if union_dtype else "cat[Record.name]",
    ).save()
    target = ln.Record(name="incoming_export_target").save()
    expected = target.name
    if union_dtype:
        label = ln.ULabel(name="incoming_export_label").save()
        expected = [target.name, label.name]
    host = ln.Record(
        name="incoming_export_host", features={feature.name: expected}
    ).save()
    referrer = ln.Record(
        name="incoming_export_referrer", features={feature.name: host.name}
    ).save()
    include = ["links_record__feature__name", "links_record__record__name"]
    if include_feature_id:
        include.append("links_record__feature_id")

    frame = ln.Record.filter(id=host.id).to_dataframe(
        features=[feature.name], include=include.copy(), is_run_input=False
    )

    actual = frame[feature.name].iloc[0]
    if union_dtype:
        assert actual == set(expected)
    else:
        assert actual == expected
    assert frame["links_record__feature__name"].iloc[0] == {feature.name}
    assert frame["links_record__record__name"].iloc[0] == {referrer.name}
    if include_feature_id:
        assert frame["links_record__feature_id"].iloc[0] == {feature.id}
    assert set(frame.columns) == {"uid", "name", feature.name, *include}
