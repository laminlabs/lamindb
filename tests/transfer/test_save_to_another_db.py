import subprocess
import sys

import lamindb as ln
import pytest


def test_cross_instance_m2m_add_still_blocked():
    # the cross-instance relation router must NOT let `.add()` link two
    # already-saved records that live on different instances; that guardrail
    # (a friendly ValueError) must still fire. Mirrors the network-dependent
    # tests/pydata test_unsaved_model_different_instance, but runs locally.
    using = f"{ln.setup.settings.user.handle}/testdb1"
    db2 = ln.DB(using)
    af = db2.Artifact.get(key="README.md")  # lives on testdb1 (see conftest)
    assert af._state.db == using

    new_label = ln.Record(name="guardrail-testlabel").save()  # on default testdb2
    try:
        with pytest.raises(ValueError) as excinfo:
            af.records.add(new_label)
        assert "Cannot label a record from instance" in str(excinfo.value)
        assert using in str(excinfo.value)
    finally:
        new_label.delete(permanent=True)


def test_lamindb_router_registered_on_fresh_import():
    code = (
        "import lamindb as ln\n"
        "from django.conf import settings\n"
        "print('LaminDBRouter' in str(settings.DATABASE_ROUTERS))\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True
    )
    assert "True" in result.stdout, result.stdout + result.stderr


def test_construct_relation_to_remote_type_before_any_save():
    using = f"{ln.setup.settings.user.handle}/testdb1"
    db2 = ln.DB(using)
    type_name = "router-remote-type"

    if db2.Record.filter(name=type_name).exists():
        db2.Record.filter(name=type_name).delete(permanent=True)
    try:
        remote_type = ln.Record(name=type_name, is_type=True).save(using=using)
        # read it back so the FK target is a record whose _state.db is the alias
        remote_type = db2.Record.get(name=type_name)
        # construction with a cross-instance FK must not raise
        child = ln.Record(name="router-remote-child", type=remote_type)
        assert child.type.uid == remote_type.uid
    finally:
        db2.Record.filter(name=type_name).delete(permanent=True)


def test_save_ulabel_to_another_db_via_model_save():
    assert ln.setup.settings.instance.name == "testdb2"

    using = f"{ln.setup.settings.user.handle}/testdb1"
    db2 = ln.DB(using)
    name = "save-using-single"
    remote_labels = db2.ULabel.filter(name=name)
    if remote_labels.exists():
        remote_labels.delete(permanent=True)

    try:
        ulabel = ln.ULabel(name=name).save(using=using)

        assert ulabel._state.db == using
        assert db2.ULabel.get(name=name)
        assert ln.ULabel.filter(name=name).count() == 0
    finally:
        db2.ULabel.filter(name=name).delete(permanent=True)


def test_save_ulabels_to_another_db_via_ln_save():
    assert ln.setup.settings.instance.name == "testdb2"

    using = f"{ln.setup.settings.user.handle}/testdb1"
    db2 = ln.DB(using)
    names = [f"save-using-bulk-{i}" for i in range(3)]
    remote_labels = db2.ULabel.filter(name__in=names)
    if remote_labels.exists():
        remote_labels.delete(permanent=True)

    try:
        ulabels = [ln.ULabel(name=name) for name in names]
        ln.save(ulabels, using=using)

        remote_matches = db2.ULabel.filter(name__in=names)
        assert remote_matches.count() == len(names)
        assert ln.ULabel.filter(name__in=names).count() == 0
    finally:
        db2.ULabel.filter(name__in=names).delete(permanent=True)


def test_save_record_with_relation_to_another_db_via_model_save():
    assert ln.setup.settings.instance.name == "testdb2"

    using = f"{ln.setup.settings.user.handle}/testdb1"
    db2 = ln.DB(using)
    type_name = "save-using-relation-type"
    rec_name = "save-using-relation-child"

    for qs in (db2.Record.filter(name=rec_name), db2.Record.filter(name=type_name)):
        if qs.exists():
            qs.delete(permanent=True)

    try:
        rec_type = ln.Record(name=type_name, is_type=True).save(using=using)
        child = ln.Record(name=rec_name, type=rec_type).save(using=using)

        assert child._state.db == using
        assert db2.Record.get(name=rec_name)
        assert ln.Record.filter(name=rec_name).count() == 0
    finally:
        for qs in (db2.Record.filter(name=rec_name), db2.Record.filter(name=type_name)):
            qs.delete(permanent=True)

def test_save_record_with_multivalued_relation_to_another_db_via_model_save():
    assert ln.setup.settings.instance.name == "testdb2"

    using = f"{ln.setup.settings.user.handle}/testdb1"
    db2 = ln.DB(using)

    type_name = "mv-intervention-type"
    gene_type_name = "mv-gene-type"
    gene_feat_name = "mv-genes"
    gene_names = [f"mv-gene-{i}" for i in range(5)]
    rec_name = "mv-child"

    def _clean():
        for qs in (
            db2.Record.filter(name=rec_name),
            db2.Record.filter(name__in=gene_names),
            db2.Record.filter(name=type_name),
            db2.Record.filter(name=gene_type_name),
            db2.Feature.filter(name=gene_feat_name),
        ):
            if qs.exists():
                qs.delete(permanent=True)

    _clean()
    try:
        intervention_t = ln.Record(name=type_name, is_type=True).save(using=using)
        gene_t = ln.Record(name=gene_type_name, is_type=True).save(using=using)
        gene_feat = ln.Feature(name=gene_feat_name, dtype=gene_t).save(using=using)
        gene_pool = [ln.Record(name=g, type=gene_t).save(using=using) for g in gene_names]

        child = ln.Record(
            name=rec_name,
            type=intervention_t,
            features={gene_feat: gene_pool[:3]},
        ).save(using=using)

        assert child._state.db == using
        assert db2.Record.get(name=rec_name)
        assert ln.Record.filter(name=rec_name).count() == 0
    finally:
        _clean()


def test_save_schema_with_features_to_another_db_via_model_save():
    assert ln.setup.settings.instance.name == "testdb2"

    using = f"{ln.setup.settings.user.handle}/testdb1"
    db2 = ln.DB(using)

    feat_names = [f"save-using-schema-feat-{i}" for i in range(3)]
    schema_name = "save-using-schema"

    def _clean():
        for qs in (
            db2.Schema.filter(name=schema_name),
            db2.Feature.filter(name__in=feat_names),
        ):
            if qs.exists():
                qs.delete(permanent=True)

    _clean()
    try:
        features = [
            ln.Feature(name=name, dtype=str).save(using=using) for name in feat_names
        ]
        schema = ln.Schema(name=schema_name, features=features).save(using=using)

        assert schema._state.db == using
        remote_schema = db2.Schema.get(name=schema_name)
        assert set(remote_schema.members.to_list("name")) == set(feat_names)
        assert ln.Schema.filter(name=schema_name).count() == 0
        assert ln.Feature.filter(name__in=feat_names).count() == 0
    finally:
        _clean()


def test_save_schema_with_remote_features_infers_instance_via_bare_save():
    assert ln.setup.settings.instance.name == "testdb2"

    using = f"{ln.setup.settings.user.handle}/testdb1"
    db2 = ln.DB(using)

    feat_names = [f"save-using-schema-infer-feat-{i}" for i in range(3)]
    schema_name = "save-using-schema-infer"

    def _clean():
        for qs in (
            db2.Schema.filter(name=schema_name),
            db2.Feature.filter(name__in=feat_names),
        ):
            if qs.exists():
                qs.delete(permanent=True)

    _clean()
    try:
        features = [
            ln.Feature(name=name, dtype=str).save(using=using) for name in feat_names
        ]
        # bare save: instance is inferred from the (remote) member features
        schema = ln.Schema(name=schema_name, features=features).save()

        assert schema._state.db == using
        remote_schema = db2.Schema.get(name=schema_name)
        assert set(remote_schema.members.to_list("name")) == set(feat_names)
        assert ln.Schema.filter(name=schema_name).count() == 0
        assert ln.Feature.filter(name__in=feat_names).count() == 0
    finally:
        _clean()


def test_save_record_with_single_valued_feature_to_another_db_via_model_save():
    assert ln.setup.settings.instance.name == "testdb2"

    using = f"{ln.setup.settings.user.handle}/testdb1"
    db2 = ln.DB(using)

    gene_type_name = "sv-gene-type"
    gene_feat_name = "sv-gene"
    gene_name = "sv-gene-value"
    rec_name = "sv-child"

    def _clean():
        for qs in (
            db2.Record.filter(name=rec_name),
            db2.Record.filter(name=gene_name),
            db2.Record.filter(name=gene_type_name),
            db2.Feature.filter(name=gene_feat_name),
        ):
            if qs.exists():
                qs.delete(permanent=True)

    _clean()
    try:
        gene_t = ln.Record(name=gene_type_name, is_type=True).save(using=using)
        gene_feat = ln.Feature(name=gene_feat_name, dtype=gene_t).save(using=using)
        gene = ln.Record(name=gene_name, type=gene_t).save(using=using)

        child = ln.Record(name=rec_name, features={gene_feat: gene}).save(using=using)

        assert child._state.db == using
        remote_child = db2.Record.get(name=rec_name)
        assert remote_child.features.get_values()[gene_feat_name] == gene_name
        assert ln.Record.filter(name=rec_name).count() == 0
    finally:
        _clean()


def test_bulk_save_records_with_multivalued_features_to_another_db():
    assert ln.setup.settings.instance.name == "testdb2"

    using = f"{ln.setup.settings.user.handle}/testdb1"
    db2 = ln.DB(using)

    gene_type_name = "bulk-mv-gene-type"
    int_type_name = "bulk-mv-int-type"
    gene_feat_name = "bulk-mv-genes"
    schema_name = "bulk-mv-schema"
    gene_names = [f"bulk-mv-gene-{i}" for i in range(5)]
    rec_names = [f"bulk-mv-child-{i}" for i in range(3)]

    def _clean():
        for qs in (
            db2.Record.filter(name__in=rec_names),
            db2.Record.filter(name__in=gene_names),
            db2.Record.filter(name=int_type_name),
            db2.Schema.filter(name=schema_name),
            db2.Feature.filter(name=gene_feat_name),
            db2.Record.filter(name=gene_type_name),
        ):
            if qs.exists():
                qs.delete(permanent=True)

    _clean()
    try:
        gene_t = ln.Record(name=gene_type_name, is_type=True).save(using=using)
        gene_feat = ln.Feature(name=gene_feat_name, dtype=gene_t).save(using=using)
        schema = ln.Schema(name=schema_name, features=[gene_feat]).save(using=using)
        int_t = ln.Record(
            name=int_type_name, is_type=True, schema=schema
        ).save(using=using)
        gene_pool = [
            ln.Record(name=g, type=gene_t).save(using=using) for g in gene_names
        ]

        records = [
            ln.Record(name=rn, type=int_t, features={gene_feat: gene_pool[:3]})
            for rn in rec_names
        ]
        ln.save(records, using=using)

        for rn in rec_names:
            remote = db2.Record.get(name=rn)
            assert set(remote.features.get_values()[gene_feat_name]) == set(
                gene_names[:3]
            )
        assert ln.Record.filter(name__in=rec_names).count() == 0
    finally:
        _clean()


def test_read_features_from_another_db_via_to_dataframe():
    # read-path gap: `.to_dataframe(include="features")` on a remote queryset must
    # resolve the record type on the queryset's instance, not the default.
    assert ln.setup.settings.instance.name == "testdb2"

    using = f"{ln.setup.settings.user.handle}/testdb1"
    db2 = ln.DB(using)

    gene_type_name = "rd-gene-type"
    int_type_name = "rd-int-type"
    gene_feat_name = "rd-genes"
    schema_name = "rd-schema"
    gene_names = [f"rd-gene-{i}" for i in range(4)]
    rec_name = "rd-child"

    def _clean():
        for qs in (
            db2.Record.filter(name=rec_name),
            db2.Record.filter(name__in=gene_names),
            db2.Record.filter(name=int_type_name),
            db2.Schema.filter(name=schema_name),
            db2.Feature.filter(name=gene_feat_name),
            db2.Record.filter(name=gene_type_name),
        ):
            if qs.exists():
                qs.delete(permanent=True)

    _clean()
    try:
        gene_t = ln.Record(name=gene_type_name, is_type=True).save(using=using)
        gene_feat = ln.Feature(name=gene_feat_name, dtype=gene_t).save(using=using)
        schema = ln.Schema(name=schema_name, features=[gene_feat]).save(using=using)
        int_t = ln.Record(
            name=int_type_name, is_type=True, schema=schema
        ).save(using=using)
        gene_pool = [
            ln.Record(name=g, type=gene_t).save(using=using) for g in gene_names
        ]
        ln.Record(
            name=rec_name, type=int_t, features={gene_feat: gene_pool[:3]}
        ).save(using=using)

        # cross-instance read with feature reassembly must not raise
        df = db2.Record.filter(name=rec_name).to_dataframe(include="features")
        assert len(df) == 1
        assert gene_feat_name in df.columns
    finally:
        _clean()
