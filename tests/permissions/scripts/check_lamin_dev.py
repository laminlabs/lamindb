import subprocess

import lamindb as ln
import pytest
from lamindb_setup.core._hub_core import select_space, select_storage


def cleanup(records):
    for record in records:
        try:
            if isinstance(record, ln.Storage):
                record.artifacts.all().delete(permanent=True)
            record.delete(permanent=True)
        except Exception as e:
            print(f"Failed deleting {record}: {e}")


assert ln.setup.settings.user.handle == "testuser1"

ln.connect("laminlabs/lamin-dev")

assert ln.setup.settings.instance.slug == "laminlabs/lamin-dev"

# check that the rename resolves correctly (it was renamed)
assert ln.Artifact.connect("laminlabs/lamin-dev1072025").db == "default"

space_name = "Our test space for CI"
space = ln.Space.get(name=space_name)

# check that we throw an error if no storage location is managed by the space
storage_loc = ln.Storage.filter(space=space).one_or_none()
assert storage_loc is None, "there should be no storage location for the space yet"
with pytest.raises(ln.errors.NoStorageLocationForSpace) as error:
    ln.track(space=space_name)  # this fails to save the env artifact

# now create the storage location in the space
storage_loc = ln.Storage("create-s3", space=space).save()
ln.track(space=space_name)
try:
    assert ln.context.space.name == space_name
    ulabel = ln.ULabel(name="My test ulabel in test space").save()
    artifact = ln.Artifact(".gitignore", key="mytest").save()

    # check that exist
    ln.ULabel.get(name="My test ulabel in test space")
    ln.Artifact.get(key="mytest")

    assert ulabel.space == space  # ulabel should end up in the restricted space
    assert artifact.space == space
    # the below check doesn't work: another worker might have associated another storage location with the space, and then the artifact ends up in that
    # assert artifact.storage == storage_loc
    # hence this check
    assert artifact.storage in ln.Storage.filter(space=space)
    assert ln.context.transform.space == space
    assert ln.context.run.space == space

    # update the space of the storage location
    space2 = ln.Space.get(name="Our test space for CI 2")
    storage_loc.space = space2
    storage_loc.save()

    response_storage = select_storage(lnid=storage_loc.uid)
    response_space = select_space(lnid=space2.uid)
    assert response_storage["space_id"] == response_space["id"]

    # connect to the instance before saving
    subprocess.run(  # noqa: S602
        "lamin connect laminlabs/lamin-dev",
        shell=True,
        check=True,
    )
    result = subprocess.run(  # noqa: S602
        "lamin save .gitignore --key mytest --space 'Our test space for CI 2'",
        shell=True,
        capture_output=True,
    )
    assert "key='mytest'" in result.stdout.decode()
    assert "storage path:" in result.stdout.decode()
    assert result.returncode == 0

finally:
    cleanup(
        (
            ulabel,
            artifact,
            ln.context.transform.latest_run,
            ln.context.transform,
            storage_loc,
        )
    )
