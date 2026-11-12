import lamindb as ln
import pandas as pd
import pytest
from django.db import transaction
from lamindb.errors import MultipleResultsFound, ValidationError
from lamindb.models.schema import KNOWN_SCHEMAS


@pytest.fixture
def schema_hash_records():
    with transaction.atomic():
        try:
            feature = ln.Feature(name="scoped_hash_feature", dtype=str).save()
            type1 = ln.Schema(name="ScopedHashType1", is_type=True).save()
            type2 = ln.Schema(name="ScopedHashType2", is_type=True).save()
            yield feature, type1, type2
        finally:
            transaction.set_rollback(True)


def test_schema_hash_lookup_omitted_and_explicit_type(schema_hash_records):
    feature, type1, type2 = schema_hash_records
    first = ln.Schema([feature], type=type1, name="first scoped schema").save()
    assert ln.Schema([feature]).uid == first.uid
    assert ln.Schema([feature], type=type1).uid == first.uid
    df = pd.DataFrame({feature.name: ["value"]})
    assert ln.Schema.from_dataframe(df).uid == first.uid

    second = ln.Schema([feature], type=type2).save()
    root = ln.Schema([feature], type=None).save()
    assert len({first.uid, second.uid, root.uid}) == 3
    assert first.hash == second.hash == root.hash
    assert root.type_id is None
    assert ln.Schema([feature], type=None).uid == root.uid
    assert ln.Schema([feature], type=type2).uid == second.uid

    # Even a matching name cannot disambiguate the initial hash lookup.
    with pytest.raises(MultipleResultsFound):
        ln.Schema([feature], name=first.name, description="must not be applied")
    first.refresh_from_db()
    assert first.type_id == type1.id
    assert first.description is None

    # Inferred feature sets retain the same unscoped lookup contract.
    with pytest.raises(MultipleResultsFound):
        ln.Schema.from_values([feature.name])
    with pytest.raises(MultipleResultsFound):
        ln.Schema.from_dataframe(df)
    inferred = ln.Schema([feature], lookup=False)
    assert inferred._state.adding
    assert inferred.hash == first.hash


def test_schema_hash_lookup_only_matches_immediate_type(schema_hash_records):
    feature, type1, _ = schema_hash_records
    subtype = ln.Schema(name="ScopedHashSubtype", type=type1, is_type=True).save()
    child = ln.Schema([feature], type=subtype).save()
    parent = ln.Schema([feature], type=type1).save()
    assert child.uid != parent.uid
    assert ln.Schema([feature], type=type1).uid == parent.uid
    assert ln.Schema([feature], type=subtype).uid == child.uid


def test_schema_hash_lookup_rejects_duplicates_within_type(schema_hash_records):
    feature, type1, _ = schema_hash_records
    ln.Schema([feature], type=type1).save()
    ln.Schema([feature], type=type1, lookup=False).save()
    with pytest.raises(MultipleResultsFound):
        ln.Schema([feature], type=type1)


def test_schema_hash_lookup_excludes_trash(schema_hash_records):
    feature, type1, _ = schema_hash_records
    first = ln.Schema([feature], type=type1).save()
    first.delete(permanent=False)
    replacement = ln.Schema([feature], type=type1).save()
    assert replacement.uid != first.uid
    assert replacement.hash == first.hash
    assert ln.Schema([feature]).uid == replacement.uid


def test_schema_hash_lookup_rejects_unsaved_type(schema_hash_records):
    feature, type1, _ = schema_hash_records
    first = ln.Schema([feature], type=type1).save()
    unsaved_type = ln.Schema(name="UnsavedScopedHashType", is_type=True)
    with pytest.raises(ValueError, match="saved"):
        ln.Schema([feature], type=unsaved_type)
    first.refresh_from_db()
    assert first.type_id == type1.id


@pytest.mark.parametrize("typed_first", [False, True])
def test_known_schema_hash_has_distinct_typed_uids(schema_hash_records, typed_first):
    _, type1, type2 = schema_hash_records
    scopes = (type1, None, type2) if typed_first else (None, type1, type2)
    schemas = [ln.Schema(itype=ln.Feature, type=scope).save() for scope in scopes]
    root = next(schema for schema in schemas if schema.type_id is None)
    assert root.uid == KNOWN_SCHEMAS[root.hash]
    assert len({schema.hash for schema in schemas}) == 1
    assert len({schema.uid for schema in schemas}) == 3
    for scope, schema in zip(scopes, schemas, strict=True):
        assert ln.Schema(itype=ln.Feature, type=scope).uid == schema.uid
    with pytest.raises(MultipleResultsFound):
        ln.Schema(itype=ln.Feature)


@pytest.mark.parametrize("occupied_by", ["typed", "trashed"])
def test_known_schema_uid_already_occupied(schema_hash_records, occupied_by):
    _, type1, _ = schema_hash_records
    original = ln.Schema(itype=ln.Feature, type=None).save()
    assert original.uid == KNOWN_SCHEMAS[original.hash]
    if occupied_by == "typed":
        original.type = type1
        original.save()
    else:
        original.delete(permanent=False)

    root = ln.Schema(itype=ln.Feature, type=None).save()
    assert root.uid != original.uid
    assert root.hash == original.hash
    assert root.type_id is None
    assert ln.Schema(itype=ln.Feature, type=None).uid == root.uid


@pytest.mark.parametrize("scope", ["root", "same_type", "different_type"])
def test_schema_hash_mutation_checks_type_scope(schema_hash_records, scope):
    feature, type1, type2 = schema_hash_records
    first_type = None if scope == "root" else type1
    second_type = type2 if scope == "different_type" else first_type
    first = ln.Schema([feature], type=first_type).save()
    second = ln.Schema([feature], type=second_type, coerce=True).save()
    original_hash = first.hash
    first.coerce = True
    if scope == "different_type":
        first.save()
        first.refresh_from_db()
        assert first.hash == second.hash
        assert first.uid != second.uid
        assert first.type_id == type1.id
        assert ln.Schema([feature], type=type1, coerce=True).uid == first.uid
    else:
        with pytest.raises(ValidationError, match="same type scope"):
            first.save()
        first.refresh_from_db()
        assert first.hash == original_hash
        assert not first.coerce


def test_schema_type_move_preserves_hash_and_identity(schema_hash_records):
    feature, type1, type2 = schema_hash_records
    first = ln.Schema([feature], type=type1).save()
    original_hash, original_uid = first.hash, first.uid
    first.type = type2
    first.save()
    first.refresh_from_db()
    assert (first.hash, first.uid) == (original_hash, original_uid)
    assert ln.Schema([feature], type=type2).uid == original_uid

    second = ln.Schema([feature], type=type1).save()
    second.type = type2
    second.save()
    # Direct moves do not silently merge distinct records.
    assert second.uid != first.uid
    with pytest.raises(MultipleResultsFound):
        ln.Schema([feature], type=type2)


def test_composite_hash_lookup_remains_structural(schema_hash_records):
    feature, type1, type2 = schema_hash_records
    child1 = ln.Schema([feature], type=type1).save()
    child2 = ln.Schema([feature], type=type2).save()
    parent1 = ln.Schema(slots={"obs": child1}, otype="AnnData", type=type1).save()
    parent2 = ln.Schema(slots={"obs": child2}, otype="AnnData", type=type2).save()
    assert parent1.hash == parent2.hash
    assert parent1.uid != parent2.uid
    assert parent1.slots["obs"].uid == child1.uid
    assert parent2.slots["obs"].uid == child2.uid

    reused = ln.Schema(slots={"obs": child2}, otype="AnnData", type=type1).save()
    assert reused.uid == parent1.uid
    assert reused.slots["obs"].uid == child1.uid
    with pytest.raises(MultipleResultsFound):
        ln.Schema(slots={"obs": child1}, otype="AnnData")
