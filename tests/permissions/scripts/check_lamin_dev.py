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
assert ulabel.space.name == space_name  # ulabel should end up in the restricted space
ulabel.delete(
    permanent=True
)  # delete silently passes in case another worker deleted the ulabel
assert (
    ln.context.transform.space.name == space_name
)  # transform and run in restricted space
assert ln.context.run.space.name == space_name  # transform and run in restricted space
ln.context.transform.delete(permanent=True)
storage_loc.delete()
