import os

import lamindb as ln
import pytest
from django.db import IntegrityError


@pytest.mark.parametrize(
    "model_class,extra_kwargs",
    [
        (ln.Record, {}),
        (ln.Feature, {"dtype": "str"}),
        (ln.Schema, {"itype": ln.Feature}),
        (ln.Project, {}),
        (ln.Reference, {}),
        (ln.ULabel, {}),
    ],
)
def test_invalid_type(model_class, extra_kwargs):
    # also see test_invalid_type_record_with_schema in test_record.py
    model_name = model_class.__name__.lower()

    no_type = model_class(name="no_type", **extra_kwargs).save()
    if model_name == "schema":
        extra_kwargs["is_type"] = True  # to avoid triggering hash look up
    with pytest.raises(ValueError) as error:
        model_class(name="WithInvalidType", type=no_type, **extra_kwargs).save()
    assert error.exconly().startswith(
        f"ValueError: You can only assign a {model_name} with `is_type=True` as `type` to another {model_name}"
    )
    # test at the database level
    if os.getenv("LAMINDB_TEST_DB_VENDOR") != "sqlite":
        no_type.is_type = True
        with pytest.raises(IntegrityError) as error:
            model_class(name="WithInvalidType", type=no_type, **extra_kwargs).save()
        assert f"{model_name}_type_is_valid_fk" in error.exconly()
    no_type.delete(permanent=True)


@pytest.mark.skipif(
    os.getenv("LAMINDB_TEST_DB_VENDOR") == "sqlite", reason="Postgres-only"
)
@pytest.mark.parametrize("model_class", [ln.Record, ln.ULabel])
def test_prevent_type_cycle(model_class):
    type_a = model_class(name="TypeA", is_type=True).save()
    type_b = model_class(name="TypeB", is_type=True).save()

    # Set A's parent to B
    type_a.type = type_b
    type_a.save()  # A → B, this is fine

    # Try to set B's parent to A (would create cycle B → A → B)
    type_b.type = type_a

    with pytest.raises(Exception) as exc_info:
        type_b.save()

    assert "cycle" in str(exc_info.value).lower()

    # Try to set type to itself
    type_a.type = type_a

    with pytest.raises(Exception) as exc_info:
        type_a.save()

    assert "cycle" in str(exc_info.value).lower()

    type_a.delete(permanent=True)
    type_b.delete(permanent=True)


@pytest.mark.parametrize("model_class", [ln.Record, ln.ULabel, ln.Project])
def test_query_sub_types_super_types_instances(model_class):
    model_name = model_class.__name__.lower()

    # Create type hierarchy
    type1 = model_class(name="Type1", is_type=True).save()
    type2 = model_class(name="Type2", is_type=True, type=type1).save()
    type3 = model_class(name="Type3", is_type=True, type=type2).save()

    # Create instances
    instance1 = model_class(name=f"{model_name}1", type=type1).save()
    instance2 = model_class(name=f"{model_name}2", type=type3).save()
    instance3 = model_class(name=f"{model_name}3", type=type3).save()

    # Get the query method dynamically
    query_method = getattr(type1, f"query_{model_name}s")

    # Children
    assert getattr(type1, model_name + "s").count() == 2  # direct instances
    assert query_method().count() == 5

    # Super types
    super_types = instance3.query_types()
    assert len(super_types) == 3
    assert super_types[0] == type3
    assert super_types[1] == type2
    assert super_types[2] == type1

    # Move type2 to trash
    type2.delete()
    assert query_method().count() == 1

    # Cleanup
    instance1.delete(permanent=True)
    instance2.delete(permanent=True)
    instance3.delete(permanent=True)
    type3.delete(permanent=True)
    type2.delete(permanent=True)
    type1.delete(permanent=True)
