import lamindb as ln
import pytest
from lamindb.models._transfer import (
    normalize_transfer_config,
    transfer_notes,
    transfer_record_feature_values,
)


def test_transfer():
    db1 = ln.DB(f"{ln.setup.settings.user.handle}/testdb1")
    artifact = db1.Artifact.get(key="README.md")
    artifact.save()
    assert artifact.key == "README.md"
    assert artifact.run is not None
    assert artifact.run.finished_at is not None
    assert artifact.run.status == "completed"
    assert artifact.run.started_at is not None
    assert ln.setup.settings.storage.root_as_str.endswith("testdb2")
    assert artifact.storage.root.endswith("testdb1")


def test_schema_transfer_ulabel_dtype():
    user_handle = ln.setup.settings.user.handle

    ln.connect("testdb1")
    perturbation_type = ln.ULabel(name="PerturbationTransferTest", is_type=True).save()
    ln.ULabel(name="DMSO", type=perturbation_type).save()
    ln.ULabel(name="IFNG", type=perturbation_type).save()
    perturbation_feature = ln.Feature(
        name="perturbation", dtype=perturbation_type
    ).save()
    schema_uid = (
        ln.Schema(
            name="transfer_schema_ulabel_perturbation",
            features=[perturbation_feature],
        )
        .save()
        .uid
    )

    ln.connect("testdb2")
    db1 = ln.DB(f"{user_handle}/testdb1")
    remote_schema = db1.Schema.get(schema_uid)
    remote_member_names = remote_schema.members.to_list("name")
    assert len(remote_member_names) > 0

    existing_local = ln.Schema.filter(uid=schema_uid).one_or_none()
    if existing_local is not None:
        existing_local.delete(permanent=True)

    existing_local_perturbation_type = ln.ULabel.filter(
        name="PerturbationTransferTest", is_type=True
    ).one_or_none()
    if existing_local_perturbation_type is not None:
        ln.ULabel.filter(type=existing_local_perturbation_type).delete(permanent=True)
        existing_local_perturbation_type.delete(permanent=True)

    assert ln.ULabel.filter(name="PerturbationTransferTest", is_type=True).count() == 0
    assert ln.ULabel.filter(type__name="PerturbationTransferTest").count() == 0

    record_only = db1.Schema.get(schema_uid).save(transfer="record")
    assert record_only.members.count() == 0

    transferred = db1.Schema.get(schema_uid).save()
    assert transferred.members.to_list("name") == remote_member_names

    perturbation = transferred.members.get(name="perturbation")
    assert perturbation.dtype_as_str.startswith("cat[ULabel[")
    assert perturbation.dtype_as_object is not None
    assert perturbation.dtype_as_object.name == "PerturbationTransferTest"

    # The dtype only needs the type. DMSO and IFNG stay on the source.
    assert ln.ULabel.filter(name="PerturbationTransferTest", is_type=True).count() == 1
    assert ln.ULabel.filter(type__name="PerturbationTransferTest").count() == 0

    before_count = transferred.links_feature.count()
    transferred_repeat = db1.Schema.get(schema_uid).save()
    assert transferred_repeat.id == transferred.id
    assert transferred_repeat.links_feature.count() == before_count


def test_ulabel_type_is_stubbed_with_the_label():
    user_handle = ln.setup.settings.user.handle
    type_name = "transfer_ci_ulabel_type"
    label_name = "transfer_ci_ulabel_value"

    ln.connect("testdb1")
    ulabel_type = ln.ULabel.filter(name=type_name, is_type=True).one_or_none()
    if ulabel_type is None:
        ulabel_type = ln.ULabel(
            name=type_name, is_type=True, description="type body"
        ).save()
    label = ln.ULabel.filter(name=label_name).one_or_none()
    if label is None:
        label = ln.ULabel(name=label_name, type=ulabel_type).save()

    ln.connect("testdb2")
    existing_type = ln.ULabel.filter(name=type_name, is_type=True).one_or_none()
    if existing_type is not None:
        ln.ULabel.filter(type=existing_type).delete(permanent=True)
        existing_type.delete(permanent=True)
    existing_label = ln.ULabel.filter(name=label_name).one_or_none()
    if existing_label is not None:
        existing_label.delete(permanent=True)

    db1 = ln.DB(f"{user_handle}/testdb1")
    transferred = db1.ULabel.get(uid=label.uid).save()
    assert transferred.type.uid == ulabel_type.uid
    stub = ln.ULabel.get(uid=ulabel_type.uid)
    assert stub.is_type is True
    assert stub.name == type_name
    assert stub.description is None


def test_record_type_parent_is_stubbed():
    user_handle = ln.setup.settings.user.handle
    parent_name = "transfer_ci_biosamples"
    child_name = "transfer_ci_rnasamples"
    feature_name = "transfer_ci_biosamples_dtype"
    schema_name = "transfer_ci_rnasamples_schema"

    ln.connect("testdb1")
    existing_schema = ln.Schema.filter(name=schema_name).one_or_none()
    if existing_schema is None:
        parent = ln.Record(
            name=parent_name, is_type=True, description="parent body"
        ).save()
        child = ln.Record(name=child_name, is_type=True, type=parent).save()
        feature = ln.Feature(name=feature_name, dtype=child).save()
        schema_uid = ln.Schema(name=schema_name, features=[feature]).save().uid
    else:
        schema_uid = existing_schema.uid
        child = ln.Record.get(name=child_name)
        parent = child.type
    parent_created_by_uid = parent.created_by.uid

    ln.connect("testdb2")
    for model, name in (
        (ln.Schema, schema_name),
        (ln.Feature, feature_name),
        (ln.Record, child_name),
        (ln.Record, parent_name),
    ):
        existing = model.filter(name=name).one_or_none()
        if existing is not None:
            existing.delete(permanent=True)
    db1 = ln.DB(f"{user_handle}/testdb1")
    db1.Schema.get(schema_uid).save()

    parent_on_target = ln.Record.get(uid=parent.uid)
    assert parent_on_target.is_type is True
    assert parent_on_target.name == parent_name
    assert parent_on_target.description is None
    assert parent_on_target.created_by.uid == parent_created_by_uid
    child_on_target = ln.Record.get(uid=child.uid)
    assert child_on_target.type.uid == parent.uid


def test_dtype_transfers_record_stub_and_schema():
    user_handle = ln.setup.settings.user.handle
    type_name = "transfer_ci_dtype_samples"
    data_name = "transfer_ci_dtype_sample_row"
    sheet_name = "transfer_ci_dtype_sheet_schema"
    column_name = "transfer_ci_dtype_sheet_column"
    feature_name = "transfer_ci_dtype_samplesheet"
    parent_name = "transfer_ci_dtype_parent_schema"

    ln.connect("testdb1")
    sample_type = ln.Record.filter(name=type_name, is_type=True).one_or_none()
    if sample_type is None:
        sample_type = ln.Record(
            name=type_name, is_type=True, description="type body"
        ).save()
    if ln.Record.filter(name=data_name, is_type=False).one_or_none() is None:
        ln.Record(name=data_name, type=sample_type).save()
    column = ln.Feature.filter(name=column_name).one_or_none()
    if column is None:
        column = ln.Feature(name=column_name, dtype=str).save()
    sheet = ln.Schema.filter(name=sheet_name).one_or_none()
    if sheet is None:
        sheet = ln.Schema(name=sheet_name, features=[column]).save()
    feature = ln.Feature.filter(name=feature_name).one_or_none()
    if feature is None:
        feature = ln.Feature(
            name=feature_name,
            dtype=sample_type,
            cat_filters={"is_type": True, "schema": sheet},
        ).save()
    parent = ln.Schema.filter(name=parent_name).one_or_none()
    if parent is None:
        parent = ln.Schema(name=parent_name, features=[feature]).save()
    parent_uid = parent.uid
    assert f"schema__uid={sheet.uid}" in feature._dtype_str

    ln.connect("testdb2")
    for model, name in (
        (ln.Schema, parent_name),
        (ln.Schema, sheet_name),
        (ln.Feature, feature_name),
        (ln.Feature, column_name),
        (ln.Record, data_name),
        (ln.Record, type_name),
    ):
        existing = model.filter(name=name).one_or_none()
        if existing is not None:
            existing.delete(permanent=True)

    db1 = ln.DB(f"{user_handle}/testdb1")
    db1.Schema.get(parent_uid).save()

    stub = ln.Record.get(uid=sample_type.uid)
    assert stub.is_type is True
    assert stub.name == type_name
    assert stub.description is None
    assert ln.Record.filter(name=data_name).one_or_none() is None
    transferred_sheet = ln.Schema.get(uid=sheet.uid)
    assert transferred_sheet.name == sheet_name
    assert transferred_sheet.members.get(name=column_name).uid == column.uid

    # A later transfer of the same feature must still bring a missing schema.
    transferred_sheet.delete(permanent=True)
    assert ln.Schema.filter(uid=sheet.uid).one_or_none() is None
    db1.Schema.get(parent_uid).save()
    assert ln.Schema.get(uid=sheet.uid).name == sheet_name


def test_feature_type_is_stubbed():
    user_handle = ln.setup.settings.user.handle
    type_name = "transfer_ci_experiment_view"
    feature_name = "transfer_ci_view_member"

    ln.connect("testdb1")
    feature_type = ln.Feature.filter(name=type_name, is_type=True).one_or_none()
    if feature_type is None:
        feature_type = ln.Feature(
            name=type_name, is_type=True, description="view body"
        ).save()
        ln.Feature(name=feature_name, dtype=str, type=feature_type).save()
    feature = ln.Feature.get(name=feature_name)
    created_by_uid = feature_type.created_by.uid

    ln.connect("testdb2")
    for name in (feature_name, type_name):
        existing = ln.Feature.filter(name=name).one_or_none()
        if existing is not None:
            existing.delete(permanent=True)
    db1 = ln.DB(f"{user_handle}/testdb1")
    db1.Feature.get(feature.uid).save()

    stub = ln.Feature.get(uid=feature_type.uid)
    assert stub.is_type is True
    assert stub.name == type_name
    assert stub.description is None
    assert stub.created_by.uid == created_by_uid
    assert ln.Feature.get(uid=feature.uid).type.uid == feature_type.uid


def test_record_transfer_links_same_name_feature_by_uid():
    user_handle = ln.setup.settings.user.handle
    feature_name = "transfer_ci_uid_assay"
    record_name = "transfer_ci_uid_record"
    sheet_name = "transfer_ci_uid_sheet"

    ln.connect("testdb1")
    feature = ln.Feature.filter(name=feature_name).one_or_none()
    if feature is None:
        feature = ln.Feature(name=feature_name, dtype=str).save()
    sheet = ln.Record.filter(name=sheet_name, is_type=True).one_or_none()
    if sheet is None:
        schema = ln.Schema(name=f"{sheet_name}_schema", features=[feature]).save()
        sheet = ln.Record(name=sheet_name, is_type=True, schema=schema).save()
    record = ln.Record.filter(name=record_name).one_or_none()
    if record is None:
        record = ln.Record(
            name=record_name, type=sheet, features={feature: "kept-on-uid"}
        ).save()

    ln.connect("testdb2")
    for model, name in (
        (ln.Record, record_name),
        (ln.Record, sheet_name),
        (ln.Schema, f"{sheet_name}_schema"),
    ):
        existing = model.filter(name=name).one_or_none()
        if existing is not None:
            existing.delete(permanent=True)
    for existing in ln.Feature.filter(name=feature_name):
        existing.delete(permanent=True)
    decoy = ln.Feature(name=feature_name, dtype=str).save()
    assert decoy.uid != feature.uid

    db1 = ln.DB(f"{user_handle}/testdb1")
    db1.Record.get(uid=sheet.uid).save(transfer="annotations")
    transferred = db1.Record.get(uid=record.uid).save(transfer="annotations")
    links = list(transferred.values_json.all())
    assert len(links) == 1
    assert links[0].feature.uid == feature.uid
    assert links[0].feature.uid != decoy.uid
    assert links[0].value == "kept-on-uid"


def test_schema_transfer_feature_uid_conflict_by_name():
    user_handle = ln.setup.settings.user.handle

    ln.connect("testdb1")
    source_feature = ln.Feature(name="tissue", dtype=str).save()
    source_feature_uid = source_feature.uid
    source_schema = ln.Schema(
        name="transfer_schema_feature_uid_conflict",
        features=[source_feature],
    ).save()
    schema_uid = source_schema.uid
    assert source_schema.members.get(name="tissue").uid == source_feature_uid

    ln.connect("testdb2")
    db1 = ln.DB(f"{user_handle}/testdb1")

    existing_local = ln.Schema.filter(uid=schema_uid).one_or_none()
    if existing_local is not None:
        existing_local.delete(permanent=True)

    existing_tissue_features = ln.Feature.filter(name="tissue")
    if existing_tissue_features.exists():
        existing_tissue_features.delete(permanent=True)

    local_tissue = ln.Feature(name="tissue", dtype=str).save()
    assert local_tissue.uid != source_feature_uid

    transferred = db1.Schema.get(schema_uid).save()
    transferred_tissue = transferred.members.get(name="tissue")
    assert transferred_tissue.uid == source_feature_uid
    assert transferred_tissue.uid != local_tissue.uid


SHEET_NAME = "transfer_ci_runs"
REC_NAME = "transfer-ci-run-1"
EMPTY_NAME = "transfer-ci-empty"
SAMPLE_TYPE = "transfer_ci_sample"
SAMPLE_NAME = "transfer-ci-S-1"
FEAT_VERSION = "package_version"
FEAT_ALIASES = "transfer_ci_aliases"
FEAT_QC = "transfer_ci_qc"
FEAT_TAGS = "transfer_ci_tags"
FEAT_USER = "transfer_ci_operator"
FEAT_SAMPLE = "transfer_ci_sample_ref"
FEAT_BIONTY = "transfer_ci_bionty_organism"
QC_TYPE = "transfer_ci_qc_type"
QC_OK = "transfer_ci_pass"
QC_FAIL = "transfer_ci_fail"
README = "transfer ci readme"


def _wipe_xferci() -> None:
    from lamindb.models.record import RecordRecord

    recs = ln.Record.filter(
        name__in=[REC_NAME, EMPTY_NAME, SAMPLE_NAME, SHEET_NAME, SAMPLE_TYPE]
    )
    ids = list(recs.values_list("id", flat=True))
    if ids:
        RecordRecord.filter(record_id__in=ids).delete()
        RecordRecord.filter(value_id__in=ids).delete()
    for name in (REC_NAME, EMPTY_NAME, SAMPLE_NAME, SHEET_NAME, SAMPLE_TYPE):
        ln.Record.filter(name=name).delete(permanent=True)
    ln.Schema.filter(name=f"{SHEET_NAME}_schema").delete(permanent=True)
    ln.Feature.filter(
        name__in=[
            FEAT_VERSION,
            FEAT_ALIASES,
            FEAT_QC,
            FEAT_TAGS,
            FEAT_USER,
            FEAT_SAMPLE,
            FEAT_BIONTY,
        ]
    ).delete(permanent=True)
    ln.ULabel.filter(name__in=[QC_OK, QC_FAIL]).delete(permanent=True)
    ln.ULabel.filter(name=QC_TYPE).delete(permanent=True)


def _source_on_testdb1() -> tuple[str, str]:
    """Idempotent fixture on testdb1: notes, scalars, lists, cats, User, nested Record."""
    ln.connect("testdb1")
    existing = ln.Record.filter(name=REC_NAME).one_or_none()
    if existing is not None:
        empty = ln.Record.filter(name=EMPTY_NAME).one()
        return existing.uid, empty.uid

    qc_type = ln.ULabel(name=QC_TYPE, is_type=True).save()
    pass_label = ln.ULabel(name=QC_OK, type=qc_type).save()
    fail_label = ln.ULabel(name=QC_FAIL, type=qc_type).save()
    sample_type = ln.Record(name=SAMPLE_TYPE, is_type=True).save()
    sample = ln.Record(name=SAMPLE_NAME, type=sample_type).save()
    user = ln.User.filter(id=ln.setup.settings.user.id).one()

    features = [
        ln.Feature(name=FEAT_VERSION, dtype=str).save(),
        ln.Feature(name=FEAT_ALIASES, dtype=list[str]).save(),
        ln.Feature(name=FEAT_QC, dtype=qc_type).save(),
        # qc_type is a ULabel record, not a type; list[...] is the dtype syntax.
        ln.Feature(name=FEAT_TAGS, dtype=list[qc_type]).save(),  # type: ignore[valid-type]
        ln.Feature(name=FEAT_USER, dtype=ln.User).save(),
        ln.Feature(name=FEAT_SAMPLE, dtype=sample_type).save(),
    ]
    schema = ln.Schema(name=f"{SHEET_NAME}_schema", features=features).save()
    sheet = ln.Record(name=SHEET_NAME, is_type=True, schema=schema).save()
    source = ln.Record(
        name=REC_NAME,
        type=sheet,
        features={
            FEAT_VERSION: "2.10.0",
            FEAT_ALIASES: ["alpha", "beta"],
            FEAT_QC: pass_label,
            FEAT_TAGS: [pass_label, fail_label],
            FEAT_USER: user,
            FEAT_SAMPLE: sample,
        },
    ).save()
    ln.models.RecordBlock(record=source, content=README, kind="readme").save()
    empty = ln.Record(name=EMPTY_NAME, type=sheet).save()
    assert source.features.get_values()[FEAT_VERSION] == "2.10.0"
    # Feature is created as str so source get_values() still works. The stored
    # dtype is then set to cat[bionty.Organism] so transferring it onto testdb2
    # (no bionty) raises ValueError from parse_dtype's ValidationError.
    bionty_feat = ln.Feature(name=FEAT_BIONTY, dtype=str).save()
    type(bionty_feat).objects.filter(pk=bionty_feat.pk).update(
        _dtype_str="cat[bionty.Organism]"
    )
    return source.uid, empty.uid


@pytest.mark.parametrize(
    "transfer,expect_notes,expect_features,error",
    [
        (None, False, False, None),
        ("sqlrecord", False, False, None),
        ("record", False, False, None),
        ("notes", True, False, None),
        ("annotations", True, True, None),
        ("nope", False, False, ValueError),
    ],
)
def test_record_transfer_features_opt_in(
    transfer, expect_notes, expect_features, error
):
    user_handle = ln.setup.settings.user.handle
    rec_uid, empty_uid = _source_on_testdb1()

    ln.connect("testdb2")
    _wipe_xferci()
    db1 = ln.DB(f"{user_handle}/testdb1")
    kwargs = {} if transfer is None else {"transfer": transfer}

    if error is not None:
        with pytest.raises(error, match="transfer should be one of"):
            db1.Record.get(uid=rec_uid).save(**kwargs)
        return

    source_for_type = db1.Record.get(uid=rec_uid)
    with pytest.raises(ValueError, match="Please transfer type"):
        source_for_type.save(**kwargs)
    sheet = source_for_type.type
    db1.Record.get(uid=sheet.uid).save(
        transfer="annotations" if expect_features else "sqlrecord"
    )

    transferred = db1.Record.get(uid=rec_uid).save(**kwargs)
    values = transferred.features.get_values()
    if expect_notes:
        assert transferred.notes == README
        assert README in transferred.describe(return_str=True)
    else:
        assert transferred.notes is None
    if not expect_features:
        assert not values
        assert not ln.Feature.filter(name=FEAT_VERSION).exists()
    if expect_features:
        assert values.get(FEAT_VERSION) == "2.10.0"
        assert set(values.get(FEAT_ALIASES) or []) == {"alpha", "beta"}
        qc = values.get(FEAT_QC)
        assert getattr(qc, "name", qc) == QC_OK
        tags = values.get(FEAT_TAGS) or []
        tag_names = {getattr(t, "name", t) for t in tags}
        assert {QC_OK, QC_FAIL} <= tag_names
        operator = values.get(FEAT_USER)
        assert getattr(operator, "handle", operator) == user_handle
        sample = values.get(FEAT_SAMPLE)
        assert getattr(sample, "name", sample) == SAMPLE_NAME
        sample = ln.Record.get(name=SAMPLE_NAME)
        assert sample.created_by.handle == user_handle
        assert sample.description is None
        assert not ln.Record.filter(name=EMPTY_NAME).exists()
        source_db = f"{user_handle}/testdb1"
        source = ln.Record.objects.using(source_db).get(uid=rec_uid)
        from lamindb.models.record import RecordJson, RecordULabel

        version = ln.Feature.objects.using(source_db).get(name=FEAT_VERSION)
        tags = ln.Feature.objects.using(source_db).get(name=FEAT_TAGS)
        qc_type = ln.ULabel.objects.using(source_db).get(name=QC_TYPE)
        RecordJson.objects.using(source_db).filter(
            record_id=source.pk, feature_id=version.id
        ).update(value="9.9.9")
        ln.Record.objects.using(source_db).filter(pk=source.pk).update(
            description="scalar-rerun"
        )
        extra = ln.ULabel.objects.using(source_db).create(
            name="transfer_ci_rerun", type=qc_type
        )
        RecordULabel.objects.using(source_db).create(
            record_id=source.pk, feature_id=tags.id, value_id=extra.id
        )
        try:
            again = db1.Record.get(uid=rec_uid).save(transfer="annotations")
            again_values = again.features.get_values()
            assert again.description == "scalar-rerun"
            assert again_values.get(FEAT_VERSION) == "9.9.9"
            again_tags = {
                getattr(label, "name", label)
                for label in (again_values.get(FEAT_TAGS) or [])
            }
            assert {QC_OK, QC_FAIL, "transfer_ci_rerun"} <= again_tags
        finally:
            RecordJson.objects.using(source_db).filter(
                record_id=source.pk, feature_id=version.id
            ).update(value="2.10.0")
            ln.Record.objects.using(source_db).filter(pk=source.pk).update(
                description=None
            )
            RecordULabel.objects.using(source_db).filter(
                record_id=source.pk, value_id=extra.id
            ).delete()
            ln.ULabel.objects.using(source_db).filter(pk=extra.pk).delete(
                permanent=True
            )
            RecordULabel.filter(value__name="transfer_ci_rerun").delete()
            ln.ULabel.filter(name="transfer_ci_rerun").delete(permanent=True)
        source_db = f"{user_handle}/testdb1"
        source_sample = ln.Record.objects.using(source_db).get(name=SAMPLE_NAME)
        ln.Record.objects.using(source_db).filter(uid=source_sample.uid).update(
            description="filled on rerun"
        )
        filled = db1.Record.get(uid=source_sample.uid).save(transfer="annotations")
        assert filled.description == "filled on rerun"
        empty = db1.Record.get(uid=empty_uid).save(transfer="annotations")
        assert not empty.features.get_values().get(FEAT_VERSION)
        db1.Record.get(uid=rec_uid).save(transfer="notes")
        assert ln.Record.get(uid=rec_uid).notes == README

        source_db = f"{user_handle}/testdb1"
        transfer_notes(transferred, transferred._state.db, None)
        source = db1.Record.get(uid=rec_uid)
        from lamindb.models.record import RecordJson

        bionty_feat = ln.Feature.objects.using(source_db).get(name=FEAT_BIONTY)
        if not source.values_json.filter(feature_id=bionty_feat.id).exists():
            RecordJson.objects.using(source_db).create(
                record_id=source.id, feature_id=bionty_feat.id, value="human"
            )
        with pytest.raises(ValueError, match="required schema module"):
            transfer_record_feature_values(
                transferred,
                source_db,
                source.pk,
                None,
                {"mapped": [], "transferred": [], "run": None},
            )
    else:
        assert values.get(FEAT_VERSION) is None


@pytest.mark.parametrize(
    "transfer,default_annotations,expected",
    [
        (None, False, "sqlrecord"),
        (None, True, "annotations"),
        ("record", False, "sqlrecord"),
        ("sqlrecord", False, "sqlrecord"),
        ("notes", False, "notes"),
        ("annotations", False, "annotations"),
    ],
)
def test_normalize_transfer_config(transfer, default_annotations, expected):
    assert (
        normalize_transfer_config(transfer, default_annotations=default_annotations)
        == expected
    )


def test_transfer_missing_space_errors():
    from lamindb.errors import NoWriteAccess
    from lamindb.models._transfer import update_fk_to_default_db

    missing = ln.Space(name="restricted-perturbations", uid="noattach1")
    missing.id = 99
    record = ln.Record(name="space gate")
    record.uid = "recSpace"
    record.space = missing

    with pytest.raises(NoWriteAccess, match="restricted-perturbations") as error:
        update_fk_to_default_db(
            record,
            "space",
            None,
            {"mapped": [], "transferred": [], "run": True},
        )
    message = str(error.value)
    target = ln.setup.settings.instance.slug
    assert (
        f"attach space 'restricted-perturbations' to the target database '{target}'"
        in message
    )
    assert "Record(uid='recSpace')" in str(error.value)


def test_annotation_transfer_requires_schema_module(connected_bionty, monkeypatch):
    import bionty as bt
    import lamindb_setup as ln_setup
    import pandas as pd

    feature = ln.Feature(
        name="organism_module_gate", dtype="cat[bionty.Organism]"
    ).save()
    record = ln.Record(name="module gate record").save()
    artifact = ln.Artifact.from_dataframe(
        pd.DataFrame({"a": [1]}), key="module-gate.parquet", description="module gate"
    ).save()
    schema = ln.Schema(name="organism module gate", itype=bt.Organism).save()
    artifact.schemas.add(schema, through_defaults={"slot": "var"})

    from lamindb.models.record import RecordJson

    RecordJson(record=record, feature=feature, value="human").save()

    instance = ln_setup.settings.instance
    modules = [module for module in instance.modules if module != "bionty"]
    monkeypatch.setattr(instance, "_schema_str", ",".join(modules))

    try:
        with pytest.raises(ValueError, match="schema module"):
            transfer_record_feature_values(
                record,
                "default",
                record.pk,
                None,
                {"mapped": [], "transferred": [], "run": True},
            )
        with pytest.raises(ValueError, match="sqlrecord"):
            artifact.features._add_from(
                artifact, transfer_logs={"mapped": [], "transferred": [], "run": True}
            )
    finally:
        artifact.schemas.clear()
        artifact.delete(permanent=True)
        schema.delete(permanent=True)
        record.delete(permanent=True)
        feature.delete(permanent=True)


def _source_user(uid: str, handle: str, name: str):
    user = ln.User(uid=uid, handle=handle, name=name)
    user._state.db = "laminlabs/lamindata"
    return user


def _transfer_logs():
    # A non-None run skips creating the transfer run in these unit tests.
    return {"mapped": [], "transferred": [], "run": True}


def test_map_user_annotation_uses_same_uid():
    from types import SimpleNamespace

    from lamindb.models._transfer import _map_user_annotation

    uid = "usrAnnot"
    handle = "annot-user"
    existing = ln.User.filter(uid=uid).one_or_none()
    if existing is not None:
        existing.delete(permanent=True)

    source = _source_user(uid, handle, "Annot User")
    logs = _transfer_logs()
    try:
        assert _map_user_annotation(source, feature=None, transfer_logs=logs) == handle
        assert ln.User.filter(uid=uid).one().handle == handle
        assert ln.User.get(uid=uid).id != ln.setup.settings.user.id
        # a second sync maps the existing registry row
        assert _map_user_annotation(source, feature=None, transfer_logs=logs) == handle
        assert ln.User.filter(uid=uid).count() == 1
        feature = SimpleNamespace(_dtype_str="cat[User.uid]")
        assert _map_user_annotation(source, feature=feature, transfer_logs=logs) == uid
    finally:
        saved = ln.User.filter(uid=uid).one_or_none()
        if saved is not None:
            saved.delete(permanent=True)


def test_map_user_annotation_asks_to_add_collaborator(monkeypatch):
    from django.db import ProgrammingError
    from lamindb.errors import NoWriteAccess
    from lamindb.models._transfer import _map_user_annotation

    uid = "usrDeny1"
    source = _source_user(uid, "denied-user", "Denied")

    def deny(self, *args, **kwargs):
        raise ProgrammingError(
            'new row violates row-level security policy for table "lamindb_user"'
        )

    monkeypatch.setattr(ln.User, "save", deny)
    with pytest.raises(NoWriteAccess, match="collaborator"):
        _map_user_annotation(source, feature=None, transfer_logs=_transfer_logs())
    assert ln.User.filter(uid=uid).one_or_none() is None


DEPTH_ROOT = "transfer_depth_root"
DEPTH_SUBTYPE = "transfer_depth_subtype"
DEPTH_DIRECT = "transfer_depth_direct"
DEPTH_GRANDCHILD = "transfer_depth_grandchild"
_DEPTH_DELETE_ORDER = (DEPTH_GRANDCHILD, DEPTH_DIRECT, DEPTH_SUBTYPE, DEPTH_ROOT)


def _depth_chain_uids() -> dict[str, str]:
    """Type, subtype, one direct data record, and one record under the subtype."""
    ln.connect("testdb1")
    root = ln.Record.filter(name=DEPTH_ROOT, is_type=True).one_or_none()
    if root is None:
        root = ln.Record(name=DEPTH_ROOT, is_type=True).save()
        subtype = ln.Record(name=DEPTH_SUBTYPE, is_type=True, type=root).save()
        ln.Record(name=DEPTH_DIRECT, type=root).save()
        ln.Record(name=DEPTH_GRANDCHILD, type=subtype).save()
    return {
        DEPTH_ROOT: ln.Record.get(name=DEPTH_ROOT, is_type=True).uid,
        DEPTH_SUBTYPE: ln.Record.get(name=DEPTH_SUBTYPE, is_type=True).uid,
        DEPTH_DIRECT: ln.Record.get(name=DEPTH_DIRECT).uid,
        DEPTH_GRANDCHILD: ln.Record.get(name=DEPTH_GRANDCHILD).uid,
    }


def _wipe_depth_chain(uids: dict[str, str]) -> None:
    for name in _DEPTH_DELETE_ORDER:
        record = ln.Record.filter(uid=uids[name]).one_or_none()
        if record is not None:
            record.delete(permanent=True)


def _transferred_depth_uids(uids: dict[str, str]) -> set[str]:
    return set(
        ln.Record.filter(uid__in=list(uids.values())).values_list("uid", flat=True)
    )


def test_query_records_depth_limit():
    uids = _depth_chain_uids()
    ln.connect("testdb1")
    root = ln.Record.get(uids[DEPTH_ROOT])
    one_hop = set(root.query_records(depth=1).values_list("name", flat=True))
    assert one_hop == {DEPTH_SUBTYPE, DEPTH_DIRECT}
    full = set(root.query_records().values_list("name", flat=True))
    assert full == {DEPTH_SUBTYPE, DEPTH_DIRECT, DEPTH_GRANDCHILD}


def test_transfer_depth_follows_type_children():
    user_handle = ln.setup.settings.user.handle
    uids = _depth_chain_uids()
    ln.connect("testdb2")
    db1 = ln.DB(f"{user_handle}/testdb1")

    _wipe_depth_chain(uids)
    db1.Record.get(uids[DEPTH_ROOT]).save(depth=0)
    assert _transferred_depth_uids(uids) == {uids[DEPTH_ROOT]}

    _wipe_depth_chain(uids)
    db1.Record.get(uids[DEPTH_ROOT]).save(depth=1)
    assert _transferred_depth_uids(uids) == {
        uids[DEPTH_ROOT],
        uids[DEPTH_SUBTYPE],
        uids[DEPTH_DIRECT],
    }

    _wipe_depth_chain(uids)
    db1.Record.get(uids[DEPTH_ROOT]).save(depth=2)
    assert _transferred_depth_uids(uids) == set(uids.values())
    _wipe_depth_chain(uids)


def _q_selects_id(query, hidden_id: int) -> bool:
    from django.db.models import Q

    if not isinstance(query, Q):
        return False
    for child in query.children:
        if isinstance(child, Q) and _q_selects_id(child, hidden_id):
            return True
        if (
            isinstance(child, tuple)
            and len(child) == 2
            and child[0] in {"id", "pk"}
            and child[1] == hidden_id
        ):
            return True
    return False


def _hide_record(monkeypatch, hidden_id: int):
    """Treat one record as missing from id lookups, the way row-level security does."""
    from lamindb.models.query_set import BasicQuerySet

    original = BasicQuerySet.filter

    def filter(self, *args, **kwargs):
        kwargs = dict(kwargs)
        if getattr(self.model, "__name__", None) == "Record":
            if kwargs.get("id") == hidden_id or kwargs.get("pk") == hidden_id:
                kwargs.pop("pk", None)
                kwargs["id"] = -1
            if "id__in" in kwargs:
                kept = [i for i in list(kwargs["id__in"]) if i != hidden_id]
                kwargs["id__in"] = kept or [-1]
            if any(_q_selects_id(arg, hidden_id) for arg in args):
                return original(self, id=-1)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(BasicQuerySet, "filter", filter)


def test_unreadable_annotation_blocks_transfer(monkeypatch):
    from uuid import uuid4

    from lamindb.models._transfer import _linked_feature_values
    from lamindb.models.record import RecordRecord
    from lamindb_setup.errors import NoReadAccess

    token = uuid4().hex[:8]
    treatment = ln.Feature(name=f"treatment-{token}", dtype="cat[Record]").save()
    combo = ln.Feature(name=f"combo-{token}", dtype="list[cat[Record]]").save()
    parent = ln.Record(name=f"parent-{token}").save()
    secret = ln.Record(name=f"secret-{token}").save()
    visible = ln.Record(name=f"visible-{token}").save()
    RecordRecord(record=parent, feature=treatment, value=secret).save()
    RecordRecord(record=parent, feature=combo, value=secret).save()
    RecordRecord(record=parent, feature=combo, value=visible).save()
    try:
        _hide_record(monkeypatch, secret.id)
        with pytest.raises(NoReadAccess, match="sqlrecord") as error:
            _linked_feature_values(parent)
        message = str(error.value)
        assert secret.name not in message
        assert f"Feature {treatment.name!r} (uid={treatment.uid})" in message
        assert f"id={secret.id}" in message
        assert f"Feature {combo.name!r} (uid={combo.uid})" in message
        assert "1 of 2 values" in message
        assert "incomplete" in message
    finally:
        RecordRecord.filter(record=parent).delete(permanent=True)
        parent.delete(permanent=True)
        secret.delete(permanent=True)
        visible.delete(permanent=True)
        treatment.delete(permanent=True)
        combo.delete(permanent=True)


def test_annotation_read_check_runs_before_save(monkeypatch):
    from lamindb.models._transfer import transfer_to_default_db
    from lamindb_setup.errors import NoReadAccess

    saved: list = []

    def blocked(record):
        raise NoReadAccess("hidden annotation")

    monkeypatch.setattr(
        "lamindb.models._transfer._save_transferred_record",
        lambda record: saved.append(record),
    )
    monkeypatch.setattr("lamindb.models._transfer._linked_feature_values", blocked)
    record = ln.Record(name="annotation gate")
    record.uid = "hidGateUID000001"
    record._state.db = "laminlabs/source"
    with pytest.raises(NoReadAccess, match="hidden annotation"):
        transfer_to_default_db(
            record,
            None,
            transfer_logs={"mapped": [], "transferred": [], "run": True},
            transfer_annotations=True,
        )
    assert saved == []


def test_clear_annotation_links_skips_unsaved_fresh_and_nonvalue_tables():
    from lamindb.models._transfer import _clear_annotation_links
    from lamindb.models.record import RecordULabel

    ln.connect("testdb2")
    unsaved = ln.Record(name="transfer_clear_unsaved")
    assert unsaved.pk is None
    _clear_annotation_links(unsaved, {})

    qc_type = ln.ULabel(name="transfer_clear_type", is_type=True).save()
    label = ln.ULabel(name="transfer_clear_label", type=qc_type).save()
    feature = ln.Feature(name="transfer_clear_feat", dtype=qc_type).save()
    record = ln.Record(name="transfer_clear_host").save()
    RecordULabel(record=record, feature=feature, value=label).save()
    transform = ln.Transform(key="transfer_clear_run.py").save()
    run = ln.Run(transform).save()
    try:
        # A row created in this sync has no links to replace.
        _clear_annotation_links(record, {"_inserted": {record.uid}})
        assert RecordULabel.filter(record=record, value=label).count() == 1
        # Run.values_artifact stores `artifact`, not `value`, so that table is skipped.
        _clear_annotation_links(run, {})
        _clear_annotation_links(record, {})
        assert RecordULabel.filter(record=record, value=label).count() == 0
    finally:
        RecordULabel.filter(record=record).delete()
        record.delete(permanent=True)
        run.delete(permanent=True)
        transform.delete(permanent=True)
        feature.delete(permanent=True)
        label.delete(permanent=True)
        qc_type.delete(permanent=True)


def test_depth_descendants_stops_and_walks():
    from lamindb.models._transfer import _depth_descendants

    uids = _depth_chain_uids()
    ln.connect("testdb1")
    root = ln.Record.get(uids[DEPTH_ROOT])
    assert _depth_descendants(root, 0) == []
    assert {record.name for record in _depth_descendants(root, 2)} == {
        DEPTH_SUBTYPE,
        DEPTH_DIRECT,
        DEPTH_GRANDCHILD,
    }

    user_handle = ln.setup.settings.user.handle
    ln.connect("testdb2")
    _wipe_depth_chain(uids)
    db1 = ln.DB(f"{user_handle}/testdb1")
    db1.Record.get(uids[DEPTH_ROOT]).save(transfer="annotations", depth=1)
    assert _transferred_depth_uids(uids) == {
        uids[DEPTH_ROOT],
        uids[DEPTH_SUBTYPE],
        uids[DEPTH_DIRECT],
    }
    _wipe_depth_chain(uids)


def test_transfer_helper_early_returns():
    from lamindb.models._transfer import (
        _cached_or_load,
        _pop_cached_linked_values,
        _put_entity,
        _read_link_rows,
        _registry_class_name,
        _remember_target,
        log_transferred_record,
        prime_annotation_transfer,
        sync_objects_from_database,
        transfer_notes,
    )

    ln.connect("testdb2")
    nameless = type("Nameless", (), {"uid": None})()
    assert _cached_or_load(nameless, {}) is None
    _remember_target(nameless, {})
    bucket: dict = {}
    _put_entity(bucket, nameless)
    _put_entity(bucket, "not-a-record")
    assert bucket == {}

    feature = ln.Feature(name="transfer_early_feat", dtype=str).save()
    record = ln.Record(name="transfer_early_host").save()
    try:
        assert _read_link_rows(feature) == ([], {})
        assert _pop_cached_linked_values({}, record) is None
        assert _pop_cached_linked_values({"_linked_values": []}, record) is None
        cached = {"_linked_values": {record.uid: []}}
        prime_annotation_transfer([feature, record], cached)
        assert cached["_linked_values"][record.uid] == []
        transfer_notes(record, record._state.db, None)
        transfer_notes(type("Plain", (), {})(), "default", 1)
        log_transferred_record(
            type("Keyed", (), {"key": "only-key", "uid": "u"})(),
            {"mapped": [], "transferred": ["x"]},
            0,
            0,
        )
        log_transferred_record(
            type("UidOnly", (), {"uid": "only-uid"})(),
            {"mapped": ["a"], "transferred": []},
            0,
            0,
        )
    finally:
        record.delete(permanent=True)
        feature.delete(permanent=True)

    user_handle = ln.setup.settings.user.handle
    source = f"{user_handle}/testdb1"
    assert _registry_class_name("ulabel") == "ULabel"
    with pytest.raises(ValueError, match="Unknown registry"):
        sync_objects_from_database("not a registry", ["uid"], source=source)
    with pytest.raises(ValueError, match="at least one uid"):
        sync_objects_from_database("record", [], source=source)
    with pytest.raises(TypeError, match="str or list"):
        sync_objects_from_database("record", None, source=source)


def test_depth_rejected_for_non_hastype_and_none():
    from lamindb.models._transfer import sync_objects_from_database

    user_handle = ln.setup.settings.user.handle
    ln.connect("testdb2")
    artifact = ln.DB(f"{user_handle}/testdb1").Artifact.get(key="README.md")
    with pytest.raises(ValueError, match="depth applies only"):
        artifact.save(depth=1)
    with pytest.raises(ValueError, match="depth must be an int >= 0"):
        artifact.save(depth=None)
    with pytest.raises(ValueError, match="depth must be an int >= 0"):
        sync_objects_from_database(
            "record", "not-a-uid", source=f"{user_handle}/testdb1", depth=None
        )
