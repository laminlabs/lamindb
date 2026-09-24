---
execute_via: python
---

# Transfer & sync across databases

This guide shows how to sync objects from a source database to your default database.

If you don't have a database, create one with the modules you need on the target.
Here we pass `bionty` because we'll transfer biological entities:

```bash
lamin init --modules bionty
```

Import `lamindb` and optionally run `ln.track()`:

```python
import lamindb as ln

ln.track()
```

Transfer works for any object type (`Artifact`, `Record`, `Transform`, `ULabel`, `Schema`, etc.).
Query the object on the source, then call `.save()` to sync it to your current default database.

## Sync an artifact

```python
db = ln.DB("laminlabs/lamindata")
# query the artifact on the source database
artifact = db.Artifact.get(key="example_datasets/mini_immuno/dataset1.h5ad")
# sync the artifact to the current database
artifact.save()
```

To transfer annotations, pass `transfer="annotations"`:

```python
# query again so that `artifact` points to the object on the source database
artifact = db.Artifact.get(key="example_datasets/mini_immuno/dataset1.h5ad")
# sync with annotations
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

The run that initiated the transfer is linked via `initiated_by_run`:

```python
artifact.run.initiated_by_run.transform
```

Upon calling `.save()` again, `lamindb` identifies that the object already exists in the target database and simply maps it:

```python
artifact = db.Artifact.get(key="example_datasets/mini_immuno/dataset1.h5ad")
artifact.save()
```

When you call `.save()` on an object queried from another database, you can pass `transfer`:

- `"sqlrecord"`: the object and its foreign keys
- `"notes"`: its associated notes
- `"annotations"`: its annotations

## Sync a record

A data record can be synced only after its type is already in the target database.
[EXP-RNA-032](https://lamin.ai/laminlabs/lamindata/record/mNDJgWFrkWQVW3ox) belongs to the [RNA-seq](https://lamin.ai/laminlabs/lamindata/record/gL3TbX2qZQmCwTAU) record frame, so transfer that frame first:

```python
rna_seq_frame = db.Record.get("gL3TbX2qZQmCwTAU")
rna_seq_frame.save(transfer="annotations")
```

Now transfer the experiment record:

```python
record = db.Record.get("mNDJgWFrkWQVW3ox")
record.save(transfer="annotations")
record.describe()
```

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

rna_seq = ln.Record.get("gL3TbX2qZQmCwTAU")
assert rna_seq.is_type
assert rna_seq._state.db == "default"

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
    "owner",
):
    assert got[key] == expected[key]

again = db.Record.get("mNDJgWFrkWQVW3ox").save(transfer="annotations")
assert again.id == record.id
```
