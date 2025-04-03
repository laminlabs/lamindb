import lamindb_setup as ln_setup

ln_setup.settings.auto_connect = False

import lamindb as ln

assert ln.setup.settings.user.handle == "testuser1"
ln.connect("laminlabs/lamin-dev")
space_name = "Our test space for CI"
ln.track(space=space_name)
assert ln.context.space.name == space_name
ulabel = ln.ULabel(name="My test ulabel in test space").save()
assert ulabel.space.name == space_name
ulabel.delete()  # delete silently passes in case another worker deleted the ulabel
ln.context.transform.delete()
