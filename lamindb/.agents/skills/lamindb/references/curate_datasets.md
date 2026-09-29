Use this when schema + dataset come from a source LaminDB instance.

```python
import lamindb as ln

SOURCE_INSTANCE = "..."  # slug like laminlabs/cellxgene
SCHEMA_UID = "..."  # a uid to be provided by the user
ARTIFACT_UID = "..."  # a uid to be provided by the user

ln.track()
db = ln.DB(SOURCE_INSTANCE)
schema = db.Schema.get(SCHEMA_UID).save()
source_artifact = db.Artifact.get(ARTIFACT_UID)
df = source_artifact.load()

# Apply curation fixes here before saving (rename columns, map typos, fill missing values, etc.)

# please note that this uses the from_dataframe() constructor, **not** Artifact()
artifact = ln.Artifact.from_dataframe(df, key=source_artifact.key, schema=schema).save()
ln.finish()
```

Workflow: write script, run it once, and if it raises `ln.errors.ValidationError`, edit the dataframe-fix block, then run again.
