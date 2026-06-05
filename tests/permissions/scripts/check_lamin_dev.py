import subprocess
from pathlib import Path
from unittest.mock import patch

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
if storage_loc is not None:
    ln.Run.filter(report__storage=storage_loc).delete(permanent=True)
    storage_loc.artifacts.all().delete(permanent=True)
    storage_loc.delete(permanent=True)

with pytest.raises(ln.errors.NoStorageLocationForSpace) as error:
    ln.track(space=space_name)  # this fails to save the env artifact
    ln.context._transform = None
    ln.context._run = None

# now create the storage location in the space
storage_loc = ln.Storage("create-s3", space=space).save()
ln.track(space=space_name)
# Initialized for failure-safe cleanup in `finally`: setup can fail before
# all records are created, so some values legitimately remain `None`.
ulabel = None
artifact = None
artifact_dir = None
artifact_storage_space = None
artifact_storage_space_file = None
try:
    assert ln.context.space.name == space_name
    ulabel = ln.ULabel(name="My test ulabel in test space").save()

    # cleanup if the artifact already exists
    artifact = ln.Artifact(".gitignore", key="mytest")
    if (
        artifact_cleanup := ln.Artifact.filter(hash=artifact.hash).one_or_none()
    ) is not None:
        artifact_cleanup.delete(permanent=True)

    # cleanup if the directory artifact already exists
    artifact_dir = ln.Artifact("./scripts", key="mytest-dir")
    if (
        artifact_cleanup := ln.Artifact.filter(hash=artifact_dir.hash).one_or_none()
    ) is not None:
        artifact_cleanup.delete(permanent=True)

    artifact = ln.Artifact(".gitignore", key="mytest").save()
    artifact_dir = ln.Artifact("./scripts", key="mytest-dir").save()

    # check that exist
    ln.ULabel.get(name="My test ulabel in test space")
    ln.Artifact.get(key="mytest")
    ln.Artifact.get(key="mytest-dir")

    assert ulabel.space == space  # ulabel should end up in the restricted space
    assert artifact.space == space
    # the below check doesn't work: another worker might have associated another storage location with the space, and then the artifact ends up in that
    # assert artifact.storage == storage_loc
    # hence this check
    assert artifact.storage in ln.Storage.filter(space=space)
    assert ln.context.transform.space == space
    assert ln.context.run.space == space

    # move the artifact to another storage location
    space_test_move = ln.Space.get(name="test-move")
    original_path = artifact.path
    artifact.space = space_test_move
    # cancel save
    with patch("builtins.input", return_value="x"):
        artifact.save()
    # save to the new storage location
    with patch("builtins.input", return_value="1"):
        artifact.save()
    assert artifact.space == space_test_move
    assert artifact.storage in ln.Storage.filter(space=space_test_move)
    assert not original_path.exists()
    assert artifact.path.as_posix().startswith(artifact.storage.root)
    assert artifact.path.exists()

    # move the directory artifact to another storage location
    assert artifact_dir.space == space
    assert artifact_dir.path.is_dir()
    assert artifact_dir.storage in ln.Storage.filter(space=space)
    original_path_dir = artifact_dir.path

    artifact_dir.space = space_test_move
    # save to the new storage location
    with patch("builtins.input", return_value="0"):
        artifact_dir.save()
    assert artifact_dir.space == space_test_move
    assert artifact_dir.storage in ln.Storage.filter(space=space_test_move)
    original_path_dir.fs.invalidate_cache()
    assert not original_path_dir.exists()
    assert artifact_dir.path.as_posix().startswith(artifact_dir.storage.root)
    assert artifact_dir.path.is_dir()

    # passing storage from a different space should set artifact.space
    # when no explicit space is provided
    key_storage_space = "mytest-storage-space.txt"
    if (
        artifact_cleanup := ln.Artifact.filter(key=key_storage_space).one_or_none()
    ) is not None:
        artifact_cleanup.delete(permanent=True)
    artifact_storage_space_file = Path(".artifact-storage-space-test.txt")
    artifact_storage_space_file.write_text("artifact storage-space test\n")
    storage_in_other_space = ln.Storage.filter(space=space_test_move).first()
    artifact_storage_space = ln.Artifact(
        artifact_storage_space_file,
        key=key_storage_space,
        storage=storage_in_other_space,
    ).save()
    assert artifact_storage_space.space == space_test_move
    assert artifact_storage_space.storage == storage_in_other_space

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
    try:
        storage_loc.run = None
        storage_loc.save()
    except:  # noqa
        pass
    if artifact_storage_space_file is not None:
        artifact_storage_space_file.unlink(missing_ok=True)
    cleanup(
        tuple(
            record
            for record in (
                ulabel,
                artifact,
                artifact_storage_space,
                artifact_dir,
                ln.context.transform.latest_run,
                ln.context.transform,
                storage_loc,
            )
            if record is not None
        )
    )
