import lamindb as ln

ln.context.version = "2"
ln.track("Ro1gl7n8YrdH0002")

assert ln.context.transform.version == "2"
