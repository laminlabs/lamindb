import lamindb as ln

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

# clean up
ulabel.delete(permanent=True)
artifact.delete(permanent=True)
ln.context.transform.delete(permanent=True)
storage_loc.delete()
