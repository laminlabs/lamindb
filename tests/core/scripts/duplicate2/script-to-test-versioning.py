import lamindb as ln

ln.context.version = "2"
ln.track("Ro1gl7n8YrdH0001")

assert ln.context.transform.version == "2"
