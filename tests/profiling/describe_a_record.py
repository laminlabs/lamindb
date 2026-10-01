import lamindb as ln

db = ln.DB("laminlabs/lamindata")
record = db.Record.get("mNDJgWFrkWQVW3ox")
record.describe()
