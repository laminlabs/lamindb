import lamindb_setup as ln_setup

ln_setup.settings.auto_connect = False

import lamindb as ln

assert ln.setup.settings.user.handle == "testuser1"

ln.connect("laminlabs/lamin-dev")

assert ln.setup.settings.instance.slug == "laminlabs/lamin-dev"

artifact = ln.Artifact.get(key="mytest")
assert artifact.space.name == "Our test space for CI"
artifact.delete()
