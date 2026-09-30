Use this when schema + dataset come from a source LaminDB instance.

```python
import lamindb as ln

ln.track()
db = ln.DB("...")  # a slug like laminlabs/cellxgene
schema = db.Schema.get("...").save()  # a schema uid
source_artifact = db.Artifact.get("...")  # an artifact uid
df = source_artifact.load()

# curate here before saving (rename columns, map typos, fill missing values, etc.)

artifact = ln.Artifact.from_dataframe(df, key=source_artifact.key, schema=schema).save()
ln.finish()
```

Workflow: write script, run it once, and if it raises `ln.errors.ValidationError`, improve your curation, then run again.
