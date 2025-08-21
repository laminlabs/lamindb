import lamindb as ln
from lamindb_setup.core._hub_core import select_space, select_storage

assert ln.setup.settings.user.handle == "testuser1"

ln.connect("laminlabs/lamin-dev")

assert ln.setup.settings.instance.slug == "laminlabs/lamin-dev"

# check that the rename resolves correctly (it was renamed)
assert ln.Artifact.using("laminlabs/lamin-dev1072025").db == "default"

# make a new storage location that's goverened by the space
space_name = "Our test space for CI"
space = ln.Space.get(name=space_name)
storage_loc = ln.Storage("create-s3", space=space).save()

ln.track(space=space_name)

assert ln.context.space.name == space_name
ulabel = ln.ULabel(name="My test ulabel in test space").save()
artifact = ln.Artifact(".gitignore", key="mytest").save()

# checks
assert ulabel.space == space  # ulabel should end up in the restricted space
assert artifact.space == space
assert artifact.storage == storage_loc
assert ln.context.transform.space == space
assert ln.context.run.space == space

# update the space of the storage location
space2 = ln.Space.get(name="Our test space for CI 2")
storage_loc.space = space2
storage_loc.save()

response_storage = select_storage(lnid=storage_loc.lnid)
response_space = select_space(lnid=space2.lnid)
assert response_storage["space_id"] == response_space["id"]

# clean up
ulabel.delete(permanent=True)
artifact.delete(permanent=True)
ln.context.transform.delete(permanent=True)
storage_loc.delete()
