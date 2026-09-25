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
