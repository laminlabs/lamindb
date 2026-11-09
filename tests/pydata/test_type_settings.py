import os
from uuid import uuid4

import lamindb as ln
import pandas as pd
import pytest
from django.db import IntegrityError, connection, transaction
from django.test.utils import CaptureQueriesContext
from lamindb.models.sqlrecord import SQLRecordSettings


@pytest.fixture
def no_name_search(monkeypatch):
    monkeypatch.setattr(ln.settings.creation, "search_names", False)


@pytest.fixture
def settings_space():
    space = ln.Space(name=f"type-settings-{uuid4().hex[:8]}").save()
    yield space
    space.delete(permanent=True)


@pytest.mark.parametrize(
    "model,member_kwargs",
    [
        (ln.Record, {}),
        (ln.Feature, {"dtype": "str"}),
        (ln.Schema, {"itype": ln.Feature, "lookup": False}),
        (ln.Project, {}),
        (ln.Reference, {}),
        (ln.ULabel, {}),
    ],
)
@pytest.mark.parametrize("policy", [False, True, "exact-space", "legacy-true"])
@pytest.mark.usefixtures("no_name_search")
def test_new_types_copy_settings(model, member_kwargs, policy, settings_space):
    parent = model(name="SettingsParent", is_type=True, space=settings_space)
    assert parent.settings.single_space is False
    if policy == "legacy-true":
        parent._aux = {"ss": "1"}
    else:
        parent.settings.single_space = (
            settings_space if policy == "exact-space" else policy
        )
    parent.emoji = "🧬"
    parent._aux["parent-only"] = {"values": [1]}

    child = model(name="SettingsChild", is_type=True, type=parent)
    sibling = model(name="SettingsSibling", is_type=True, type=parent)
    grandchild = model(name="SettingsGrandchild", is_type=True, type=child)
    member = model(name="SettingsMember", type=parent, **member_kwargs)

    expected = settings_space if policy == "exact-space" else policy is not False
    assert child.settings.single_space == expected
    assert grandchild.settings.single_space == expected
    assert (
        child._aux
        == sibling._aux
        == {key: value for key, value in parent._aux.items() if key != "parent-only"}
    )
    assert member._aux == {"ei": "🧬"}
    assert child._aux is not parent._aux
    assert child._aux is not sibling._aux

    child.settings.single_space = False
    assert child._aux == {"ei": "🧬"}
    assert sibling.settings.single_space == expected
    parent.settings.single_space = False
    assert sibling.settings.single_space == expected


@pytest.mark.parametrize("skip_validation", [False, True])
@pytest.mark.parametrize("child_aux", [None, {}, {"ei": "🧪", "child-only": 2}])
@pytest.mark.usefixtures("no_name_search")
def test_copy_settings_preserves_child_metadata(skip_validation, child_aux):
    parent = ln.Record(name="SettingsParent", is_type=True)
    parent.settings.single_space = True
    parent.emoji = "🧬"
    child = ln.Record(
        name="SettingsChild",
        type=parent,
        is_type=True,
        _aux=child_aux,
        _skip_validation=skip_validation,
    )
    assert child.settings.single_space is True
    assert child.emoji == (child_aux or {}).get("ei", "🧬")
    if child_aux is not None:
        assert "ss" not in child_aux
    if child_aux:
        assert child._aux["child-only"] == 2


@pytest.mark.usefixtures("no_name_search")
def test_explicit_setting_and_future_mutable_settings(monkeypatch, settings_space):
    # Declared settings can grow without changing the constructor/helper.
    monkeypatch.setattr(
        SQLRecordSettings, "_aux_fields", {"single_space": "ss", "future": "future"}
    )
    parent = ln.Record(name="SettingsParent", is_type=True)
    parent.settings.single_space = True
    parent._aux["future"] = {"values": [1]}
    child_aux = {"ss": settings_space.uid, "child-only": {"values": [2]}}
    child = ln.Record(name="SettingsChild", type=parent, is_type=True, _aux=child_aux)
    sibling = ln.Record(name="SettingsSibling", type=parent, is_type=True)
    assert child.settings.single_space == settings_space
    child._aux["future"]["values"].append(3)
    child._aux["child-only"]["values"].append(4)
    assert parent._aux["future"] == sibling._aux["future"] == {"values": [1]}
    assert child_aux == {"ss": settings_space.uid, "child-only": {"values": [2]}}
    explicit = ln.Record(
        name="SettingsExplicitDefault",
        type=parent,
        is_type=True,
        _aux={"future": False},
    )
    assert explicit._aux["future"] is False


def test_settings_copy_does_not_query_or_refresh_parent(settings_space, monkeypatch):
    parent = ln.Record(name="SettingsSnapshot", is_type=True).save()
    try:
        # An unsaved edit is the source, even if it names a missing Space UID.
        parent._aux = {"ss": "missing-space-uid"}
        branch = ln.Branch.get(1)
        with CaptureQueriesContext(connection) as queries:
            child = ln.Record(
                name="SettingsSnapshotChild",
                type=parent,
                is_type=True,
                space=settings_space,
                branch=branch,
                _search_names=False,
            )
        assert len(queries) == 0
        assert child._aux == {"ss": "missing-space-uid"}
        assert ln.Record.get(parent.id).settings.single_space is False

        # Multi-instance construction also copies without resolving Space UIDs.
        monkeypatch.setenv("LAMINDB_MULTI_INSTANCE", "true")
        with CaptureQueriesContext(connection) as queries:
            other = ln.Record(
                name="SettingsMultiInstance",
                type=parent,
                is_type=True,
                _search_names=False,
            )
        assert len(queries) == 0
        assert other._aux == child._aux
    finally:
        parent.delete(permanent=True)


def test_settings_are_only_copied_at_construction(settings_space):
    parent = ln.Record(name="SettingsLifecycleParent", is_type=True).save()
    existing = ln.Record(name="SettingsExisting", is_type=True, type=parent).save()
    child = None
    try:
        parent.settings.single_space = True
        parent.save()
        existing.refresh_from_db()
        assert existing.settings.single_space is False
        assert ln.Record.get(existing.id).settings.single_space is False
        reused = ln.Record(name=existing.name, type=parent, is_type=True)
        assert reused.id == existing.id
        assert reused.settings.single_space is False
        reused.save()

        child = ln.Record(name="SettingsOverride", type=parent, is_type=True)
        assert child.settings.single_space is True
        child.settings.single_space = False
        child.save()
        child.refresh_from_db()
        assert child.settings.single_space is False
        # Do not traverse to the grandparent when the immediate policy is off.
        grandchild = ln.Record(name="SettingsUnconstrained", type=child, is_type=True)
        assert grandchild.settings.single_space is False

        late_type = ln.Record(name="SettingsLateType", is_type=True)
        late_type.type = parent
        assert late_type.settings.single_space is False

        # An explicit True remains relative to the child, not a fixed parent UID.
        relative = ln.Record(
            name="SettingsRelative", type=parent, is_type=True, space=settings_space
        )
        assert relative._aux["ss"] == 1
        assert (
            ln.Record(
                name="SettingsRelativeMember", type=relative, _search_names=False
            ).space
            == settings_space
        )
    finally:
        if child is not None and child.id is not None:
            child.delete(permanent=True)
        existing.delete(permanent=True)
        parent.delete(permanent=True)


def test_schema_implied_type_and_dataframe_destination(settings_space):
    feature = ln.Feature(name="settings_value", dtype="int").save()
    schema = ln.Schema([feature]).save()
    parent = ln.Record(name="SettingsSheetParent", is_type=True).save()
    sheet = None
    rows = []
    try:
        parent.settings.single_space = settings_space
        parent.save()
        sheet = ln.Record(name="SettingsSheet", type=parent, schema=schema)
        assert sheet.is_type is True
        assert sheet.settings.single_space == settings_space
        assert sheet.space == settings_space
        sheet.save()
        rows = ln.Record.from_dataframe(
            pd.DataFrame({"settings_value": [1, 2]}), type=sheet
        ).save()
        assert len(rows) == 2
        assert all(row.space == settings_space for row in rows)
        assert all("ss" not in (row._aux or {}) for row in rows)
    finally:
        for row in rows:
            row.delete(permanent=True)
        if sheet is not None and sheet.id is not None:
            sheet.delete(permanent=True)
        parent.delete(permanent=True)
        schema.delete(permanent=True)
        feature.delete(permanent=True)


def test_manager_and_bulk_creation_copy_settings():
    parent = ln.Record(name="SettingsBatchParent", is_type=True).save()
    children = []
    try:
        parent.settings.single_space = True
        parent.save()
        child, created = ln.Record.objects.get_or_create(
            name="SettingsManagerChild", type=parent, is_type=True
        )
        children.append(child)
        assert created
        assert child.settings.single_space is True
        child.settings.single_space = False
        child.save()
        reused, created = ln.Record.objects.get_or_create(
            name=child.name, type=parent, is_type=True
        )
        assert not created
        assert reused.settings.single_space is False

        batch = [
            ln.Record(name=f"SettingsBatchChild{i}", type=parent, is_type=True)
            for i in range(2)
        ]
        ln.save(batch)
        children.extend(batch)
        assert all(
            ln.Record.get(record.id).settings.single_space is True for record in batch
        )
    finally:
        for child in reversed(children):
            child.delete(permanent=True)
        parent.delete(permanent=True)


@pytest.mark.skipif(
    os.getenv("LAMINDB_TEST_DB_VENDOR") == "sqlite", reason="Postgres-only"
)
def test_settings_copy_respects_parent_visibility(settings_space):
    # This suite uses a disposable admin database. Roll back the role, policies,
    # and records together so later tests retain their normal database access.
    with transaction.atomic():
        try:
            parent_space = ln.Space(name="settings-read-only-parent").save()
            parent = ln.Record(
                name="SettingsReadableParent", is_type=True, space=parent_space
            )
            parent.settings.single_space = settings_space
            parent.save()
            role = connection.ops.quote_name(f"settings_reader_{uuid4().hex[:8]}")
            with connection.cursor() as cursor:
                cursor.execute("SET CONSTRAINTS ALL IMMEDIATE")
                cursor.execute(f"CREATE ROLE {role} NOLOGIN")
                cursor.execute(f"GRANT USAGE ON SCHEMA public TO {role}")
                cursor.execute(f"GRANT SELECT ON ALL TABLES IN SCHEMA public TO {role}")
                cursor.execute(
                    f"GRANT USAGE ON ALL SEQUENCES IN SCHEMA public TO {role}"
                )
                cursor.execute(f"GRANT INSERT ON lamindb_record TO {role}")
                cursor.execute("ALTER TABLE lamindb_record ENABLE ROW LEVEL SECURITY")
                cursor.execute(
                    f"CREATE POLICY settings_read ON lamindb_record FOR SELECT TO {role} "
                    f"USING (space_id IN ({parent_space.id}, {settings_space.id}))"
                )
                cursor.execute(
                    f"CREATE POLICY settings_insert ON lamindb_record FOR INSERT TO {role} "
                    f"WITH CHECK (space_id = {settings_space.id})"
                )
                cursor.execute(f"SET LOCAL ROLE {role}")

            # Reading the parent is enough; writing to its space is unnecessary.
            loaded_parent = ln.Record.get(parent.id)
            child = ln.Record(
                name="SettingsReadableChild", type=loaded_parent, is_type=True
            )
            assert child.settings.single_space == settings_space
            child.save()
            assert ln.Record.get(child.id)._aux["ss"] == settings_space.uid

            with connection.cursor() as cursor:
                cursor.execute("RESET ROLE")
                cursor.execute(
                    f"ALTER POLICY settings_read ON lamindb_record "
                    f"USING (space_id = {settings_space.id})"
                )
                cursor.execute(f"SET LOCAL ROLE {role}")

            with pytest.raises(ln.Record.DoesNotExist):
                ln.Record.get(parent.id)

            # No fresh permission query during copying, even after revocation.
            blocked = ln.Record(
                name="SettingsHiddenParent", type=loaded_parent, is_type=True
            )
            assert blocked._aux["ss"] == settings_space.uid
            with (
                pytest.raises(IntegrityError, match="record_type_is_valid_fk"),
                transaction.atomic(),
            ):
                blocked.save()

            # A low-level insert (as used by REST) also cannot use the hidden parent.
            with (
                pytest.raises(IntegrityError, match="record_type_is_valid_fk"),
                transaction.atomic(),
            ):
                with connection.cursor() as cursor:
                    cursor.execute(
                        "INSERT INTO lamindb_record (uid, name, is_type, type_id, space_id, created_by_id) "
                        "VALUES (%s, %s, TRUE, %s, %s, %s)",
                        (
                            blocked.uid,
                            blocked.name,
                            parent.id,
                            settings_space.id,
                            child.created_by_id,
                        ),
                    )

            # The child's copied policy works without reading its hidden ancestor.
            grandchild = ln.Record(
                name="SettingsVisibleGrandchild", type=child, is_type=True
            ).save()
            assert grandchild.space == settings_space
        finally:
            transaction.set_rollback(True)
