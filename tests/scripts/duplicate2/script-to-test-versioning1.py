import lamindb as ln

ln.context.uid = "Ro1gl7n8YrdH0001"
ln.context.version = "2"

ln.context.track()

assert ln.context.transform.version == "2"
