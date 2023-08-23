import pytest
from django.db.models.deletion import ProtectedError

import lamindb as ln


def test_is_new_version_of_versioned_transform():
    # attempt to create a transform with an invalid version
    with pytest.raises(ValueError) as error:
        transform = ln.Transform(name="My transform", version=0)
    assert (
        error.exconly()
        == "ValueError: `version` parameter must be `None` or `str`, e.g., '0.1', '1',"
        " '2', etc."
    )
    # with pytest.raises(ValueError) as error:
    #     transform = ln.Transform(name="My transform", version="0")
    # assert (
    #     error.exconly()
    #     == "ValueError: Please choose a version != '0', as it could be interpreted as `None`"  # noqa
    # )

    # create a versioned transform
    transform = ln.Transform(name="My transform", version="1")
    assert transform.version == "1"

    transform.save()

    # create new transform from old transform
    transform_v2 = ln.Transform(name="My 2nd transform", is_new_version_of=transform)
    assert transform_v2.id[:12] == transform.id[:12]  # stem_id
    assert transform.version == "1"
    assert (
        transform.initial_version_id is None
    )  # initial transform has initial_version_id None
    assert transform_v2.initial_version_id == transform.id
    assert transform_v2.version == "2"

    transform_v2.save()

    # create new transform from newly versioned transform
    transform_v3 = ln.Transform(name="My transform", is_new_version_of=transform_v2)
    assert transform_v3.id[:12] == transform.id[:12]  # stem_id
    assert transform_v3.initial_version_id == transform.id
    assert transform_v3.version == "3"

    # test that reference transform cannot be deleted
    with pytest.raises(ProtectedError):
        transform.delete()
    transform_v2.delete()
    transform_v3.delete()
    transform.delete()


def test_is_new_version_of_unversioned_transform():
    # unversioned transform
    transform = ln.Transform(name="My transform")
    assert transform.initial_version_id is None
    assert transform.version is None

    # what happens if we don't save the old transform?
    # add a test for it!
    transform.save()

    # create new transform from old transform
    new_transform = ln.Transform(name="My new transform", is_new_version_of=transform)
    assert new_transform.id[:12] == transform.id[:12]  # stem_id
    assert transform.version == "1"
    assert transform.initial_version is None
    assert new_transform.initial_version_id == transform.id
    assert new_transform.version == "2"

    transform.delete()
