from unittest.mock import patch

import numpy as np
import pytest

import lamindb as ln
from lamindb.models._feature_manager import FeatureManager
from lamindb.models.sqlrecord import (
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

    perturbation_labels = sorted(
        perturbation.dtype_as_object.ulabels.values_list("name", flat=True)
    )
    assert perturbation_labels == ["DMSO", "IFNG"]

    before_count = transferred.links_feature.count()
    transferred_repeat = db1.Schema.get(schema_uid).save()
    assert transferred_repeat.id == transferred.id
    assert transferred_repeat.links_feature.count() == before_count


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
FEAT_SKIP = "transfer_ci_skipped_module"
ORG_NAME = "transfer-ci-organism"
ORG_TYPE = "transfer_ci_organism_type"
ORG_SCHEMA = "transfer_ci_organism_schema"
QC_TYPE = "transfer_ci_qc_type"
QC_PASS = "transfer_ci_pass"
QC_FAIL = "transfer_ci_fail"
README = "transfer ci readme"


def _wipe_xferci() -> None:
    from lamindb.models.record import RecordRecord

    recs = ln.Record.filter(
        name__in=[
            REC_NAME,
            EMPTY_NAME,
            SAMPLE_NAME,
            SHEET_NAME,
            SAMPLE_TYPE,
            ORG_NAME,
            ORG_TYPE,
        ]
    )
    ids = list(recs.values_list("id", flat=True))
    if ids:
        RecordRecord.filter(record_id__in=ids).delete()
        RecordRecord.filter(value_id__in=ids).delete()
    for name in (
        REC_NAME,
        EMPTY_NAME,
        SAMPLE_NAME,
        SHEET_NAME,
        SAMPLE_TYPE,
        ORG_NAME,
        ORG_TYPE,
    ):
        ln.Record.filter(name=name).delete(permanent=True)
    ln.Schema.filter(name=f"{SHEET_NAME}_schema").delete(permanent=True)
    ln.Schema.filter(name=ORG_SCHEMA).delete(permanent=True)
    ln.Feature.filter(
        name__in=[
            FEAT_VERSION,
            FEAT_ALIASES,
            FEAT_QC,
            FEAT_TAGS,
            FEAT_USER,
            FEAT_SAMPLE,
            FEAT_SKIP,
        ]
    ).delete(permanent=True)
    ln.ULabel.filter(name__in=[QC_PASS, QC_FAIL]).delete(permanent=True)
    ln.ULabel.filter(name=QC_TYPE).delete(permanent=True)


def _source_on_testdb1() -> tuple[str, str, str]:
    """Idempotent fixture on testdb1: notes, scalars, lists, cats, User, nested Record."""
    ln.connect("testdb1")
    existing = ln.Record.filter(name=REC_NAME).one_or_none()
    if existing is not None:
        empty = ln.Record.filter(name=EMPTY_NAME).one()
        org_rec = ln.Record.filter(name=ORG_NAME).one()
        return existing.uid, empty.uid, org_rec.uid

    qc_type = ln.ULabel(name=QC_TYPE, is_type=True).save()
    pass_label = ln.ULabel(name=QC_PASS, type=qc_type).save()
    fail_label = ln.ULabel(name=QC_FAIL, type=qc_type).save()
    sample_type = ln.Record(name=SAMPLE_TYPE, is_type=True).save()
    sample = ln.Record(name=SAMPLE_NAME, type=sample_type).save()
    user = ln.User.filter(id=ln.setup.settings.user.id).one()

    features = [
        ln.Feature(name=FEAT_VERSION, dtype=str).save(),
        ln.Feature(name=FEAT_ALIASES, dtype=list[str]).save(),
        ln.Feature(name=FEAT_QC, dtype=qc_type).save(),
        ln.Feature(name=FEAT_TAGS, dtype=list[qc_type]).save(),
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

    from django.apps import apps

    Organism = apps.get_model("bionty", "Organism")
    organism = Organism.objects.filter(name="human").first()
    if organism is None:
        organism = Organism(name="human")
        organism.save()
    org_feat = ln.Feature(name=FEAT_SKIP, dtype=Organism).save()
    org_schema = ln.Schema(name=ORG_SCHEMA, features=[org_feat]).save()
    org_type = ln.Record(name=ORG_TYPE, is_type=True, schema=org_schema).save()
    org_rec = ln.Record(
        name=ORG_NAME,
        type=org_type,
        features={FEAT_SKIP: organism},
    ).save()
    return source.uid, empty.uid, org_rec.uid


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
def test_record_transfer_features_opt_in(transfer, expect_notes, expect_features, error):
    user_handle = ln.setup.settings.user.handle
    rec_uid, empty_uid, org_uid = _source_on_testdb1()

    ln.connect("testdb2")
    _wipe_xferci()
    db1 = ln.DB(f"{user_handle}/testdb1")
    kwargs = {} if transfer is None else {"transfer": transfer}

    if error is not None:
        with pytest.raises(error, match="transfer should be one of"):
            db1.Record.get(uid=rec_uid).save(**kwargs)
        return

    transferred = db1.Record.get(uid=rec_uid).save(**kwargs)
    values = transferred.features.get_values()
    if expect_notes:
        assert transferred.notes == README
    else:
        assert transferred.notes is None
    if expect_features:
        assert values.get(FEAT_VERSION) == "2.10.0"
        assert set(values.get(FEAT_ALIASES) or []) == {"alpha", "beta"}
        qc = values.get(FEAT_QC)
        assert getattr(qc, "name", qc) == QC_PASS
        tags = values.get(FEAT_TAGS) or []
        tag_names = {getattr(t, "name", t) for t in tags}
        assert {QC_PASS, QC_FAIL} <= tag_names
        operator = values.get(FEAT_USER)
        assert getattr(operator, "handle", operator) == user_handle
        sample = values.get(FEAT_SAMPLE)
        assert getattr(sample, "name", sample) == SAMPLE_NAME
        empty = db1.Record.get(uid=empty_uid).save(transfer="annotations")
        assert not empty.features.get_values().get(FEAT_VERSION)
        db1.Record.get(uid=rec_uid).save(transfer="notes")
        assert ln.Record.get(uid=rec_uid).notes == README
        with pytest.raises(ValueError, match="required schema module"):
            db1.Record.get(uid=org_uid).save(transfer="annotations")

        source_db = f"{user_handle}/testdb1"
        transfer_notes(transferred, transferred._state.db, None)
        transfer_record_feature_values(transferred, source_db, None, None, {})
        source = db1.Record.get(uid=rec_uid)
        sample_obj = ln.Record.objects.using(source_db).get(name=SAMPLE_NAME)
        qc_obj = ln.ULabel.objects.using(source_db).get(name=QC_PASS)
        user = ln.User.objects.using(source_db).get(handle=user_handle)
        orig = FeatureManager.get_values

        def _values_for_prepare(self, *args, **kwargs):
            got_values = orig(self, *args, **kwargs)
            if getattr(self._host, "uid", None) != rec_uid:
                return got_values
            if self._host._state.db != source_db:
                return got_values
            got_values[FEAT_ALIASES] = np.array(["alpha", "beta"])
            got_values[FEAT_SAMPLE] = sample_obj
            got_values[FEAT_QC] = qc_obj
            if getattr(user, "name", None):
                got_values[FEAT_USER] = user.name
            return got_values

        with patch.object(FeatureManager, "get_values", _values_for_prepare):
            transfer_record_feature_values(
                transferred,
                source_db,
                source.pk,
                None,
                {"mapped": [], "transferred": [], "run": None},
            )
        again = ln.Record.get(uid=rec_uid).features.get_values()
        assert set(again.get(FEAT_ALIASES) or []) == {"alpha", "beta"}
        assert getattr(again.get(FEAT_SAMPLE), "name", again.get(FEAT_SAMPLE)) == SAMPLE_NAME
        assert getattr(again.get(FEAT_QC), "name", again.get(FEAT_QC)) == QC_PASS
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
