import lamindb as ln
import pytest
from lamindb.errors import ValidationError


def test_case1_remove_values_then_add_values():
    """
    test removing and adding values to a dataframe
    test for the github issue 793, where when we were adding values it led to adding
    values to the present one, instead of first adding and then removing
    :return:
    """

    ln.Feature(name="a", dtype=int).save()
    ln.Feature(name="b", dtype=int).save()
    artifact = ln.Artifact(
        "/Users/ishitajain/PycharmProjects/lamindb_add_vals/lamindb/tests/fixtures/df",
        key="test_scalar_cardinality_case1",
    ).save()
    artifact.features.add_values({"a": 1, "b": 2})
    assert artifact.features.get_values() == {"a": 1, "b": 2}
    with pytest.raises(ValidationError):
        artifact.features.add_values({"a": 4, "b": 5})
    # Values unchanged after failed add
    assert artifact.features.get_values() == {"a": 1, "b": 2}


def test_case2_remove_values():
    """
    testing features remove values from an artifact
    :return:
    """

    ln.Feature(name="case_1", dtype=int).save()
    ln.Feature(name="case_2", dtype=int).save()
    artifact_3 = ln.Artifact(
        "/Users/ishitajain/PycharmProjects/lamindb_add_vals/lamindb/tests/fixtures/df_3.tsv",
        key="test_scalar_cardinality_case4.tsv",
    ).save()
    artifact_3.features.add_values({"case_1": 1, "case_2": 2})
    assert artifact_3.features.get_values() == {"case_1": 1, "case_2": 2}
    artifact_3.features.remove_values()
    artifact_3.features.add_values({"case_1": 3, "case_2": 4})
    assert artifact_3.features.get_values() == {"case_1": 3, "case_2": 4}


def test_case3_set_values():
    """Case 3: set_values replaces all external feature values.

    set_values does _remove_values() then _add_values(); final state must
    be exactly the new values (a=3, b=4) with no leftovers.
    """
    ln.Feature(name="case3_a", dtype=int).save()
    ln.Feature(name="case3_b", dtype=int).save()
    artifact_1 = ln.Artifact(
        "/Users/ishitajain/PycharmProjects/lamindb_add_vals/lamindb/tests/fixtures/df_1",
        key="test_scalar_cardinality_case2",
    ).save()

    artifact_1.features.add_values({"case3_a": 1, "case3_b": 2})
    artifact_1.features.set_values({"case3_a": 3, "case3_b": 4})

    assert artifact_1.features.get_values() == {"case3_a": 3, "case3_b": 4}
