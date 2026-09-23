---
execute_via: python
---

# Transfer & sync across databases

This guide shows how to sync objects from a source database to your default database.

If you don't have a database, create one with the modules you need on the target.
Here we pass `bionty` because we'll transfer ontology labels:

```bash
lamin init --modules bionty
```

Import `lamindb` and optionally run `ln.track()`:

```python
import lamindb as ln

ln.track()
```

Syncing works for any object type (`Artifact`, `Record`, `Transform`, `ULabel`, `Schema`, etc.).
Query the object on the source, then call `.save()` to copy it into the default database.

## Sync an artifact

```python
db = ln.DB("laminlabs/lamindata")
# query the artifact on the source database
artifact = db.Artifact.get(key="example_datasets/mini_immuno/dataset1.h5ad")
# sync the artifact to the current database
artifact.save()
```

The default (`transfer="sqlrecord"`) copies the metadata row and its foreign keys.
To also copy feature & label annotations, pass `transfer="annotations"`:

```python
# query again so that `artifact` holds the object on the source database
artifact = db.Artifact.get(key="example_datasets/mini_immuno/dataset1.h5ad")
# sync the artifact to the current database, including transfer of annotations where necessary
artifact.save(transfer="annotations")
```

The artifact now has all feature & label annotations:

```python
artifact.describe()
```

The sync is zero-copy: the data itself remains in the original storage location.

```python
artifact.path
```

Data lineage indicates the source database of the sync:

```python
artifact.view_lineage()
```

The run that initiated the sync is linked via `initiated_by_run`:

```python
artifact.run.initiated_by_run.transform
```

Upon calling `.save()` again, `lamindb` identifies that the object already exists in the target database and simply maps it:

```python
artifact = db.Artifact.get(key="example_datasets/mini_immuno/dataset1.h5ad")
artifact.save()
```

**What the `transfer` argument controls.**

`.save()` on an object queried from another database accepts `transfer`:

- `"sqlrecord"`: the row and its foreign keys
- `"notes"`: also the latest readme
- `"annotations"`: also feature values and M2M labels

{class}`~lamindb.Schema` is the exception: omitting `transfer` defaults to `"annotations"` so schema members come along.

`created_by` is always remapped to the user who runs the transfer. A {class}`~lamindb.User` feature value is remapped the same way.

If the target is missing a schema module (for example you transfer a `bionty.Organism` value but did not run `lamin init --modules bionty`), that feature is skipped with a warning.

## Sync a record

The same modes apply to records. This experiment on `laminlabs/lamindata` has scalar, {class}`~lamindb.User`, {class}`~lamindb.Project`, and bionty features: [EXP-RNA-032](https://lamin.ai/laminlabs/lamindata/record/mNDJgWFrkWQVW3ox).

```python
record = db.Record.get("mNDJgWFrkWQVW3ox")
record.save()
record.features.get_values()
```

```python tags=["hide-cell"]
assert record.uid == "mNDJgWFrkWQVW3ox"
assert record._state.db == "default"
assert record.type is not None
assert record.features.get_values() == {"name": "EXP-RNA-032"}
```

The row and its type are now local. `name` is stored on the row, so it already appears; the other feature values do not. Pass `transfer="annotations"` to copy them:

```python
# query again so that `record` holds the object on the source database
record = db.Record.get("mNDJgWFrkWQVW3ox")
record.save(transfer="annotations")
```

```python
record.describe()
```

`owner` now points at the current user. `organism`, `assay`, and `project` keep their identity.

```{dropdown} How do I know if an object is in the default database or elsewhere?

Every `SQLRecord` object has an attribute `._state.db` which can take the following values:

- `None`: the object has not yet been saved to any database
- `"default"`: the object is saved on the default database instance
- `"account/name"`: the object is saved on a non-default database instance referenced by `account/name` (e.g., `laminlabs/lamindata`)

```

```python tags=["hide-cell"]
assert artifact.transform.description == "Transfer from `laminlabs/lamindata`"
assert artifact.transform.key == "__lamindb_transfer__/4XIuR0tvaiXM"
assert artifact.transform.uid == "4XIuR0tvaiXM0000"
assert artifact.run.initiated_by_run.transform.description.startswith("Transfer & sync")
assert artifact.features.slots
for schema in artifact.features.slots.values():
    _ = schema.index

source = db.Record.get("mNDJgWFrkWQVW3ox")
expected = source.features.get_values()
got = record.features.get_values()
assert record._state.db == "default"
assert set(got) == set(expected)
for key in (
    "date_of_experiment",
    "organism",
    "assay",
    "project",
    "n_samples",
    "notes",
    "name",
    "description",
):
    assert got[key] == expected[key]
assert got["owner"] == ln.setup.settings.user.handle

again = db.Record.get("mNDJgWFrkWQVW3ox").save(transfer="annotations")
assert again.id == record.id
```
