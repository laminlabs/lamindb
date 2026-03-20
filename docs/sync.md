---
execute_via: python
---

# Sync data across databases

This guide shows how to sync objects from a source database to your default database.

We need a target database:

```python
!lamin init --storage ./test-sync --modules bionty
```

Import `lamindb` and optionally run `ln.track()`:

```python
import lamindb as ln

ln.track()
```

Syncing works for any object type (`Artifact`, `Record`, `Transform`, `ULabel`, etc.). Let's sync an artifact to our current default database:

```python
db = ln.DB("laminlabs/lamindata")
artifact = db.Artifact.get(key="example_datasets/mini_immuno/dataset1.h5ad")  # query the artifact
artifact.save()  # sync the artifact to the current database
```

If you also want to sync feature & label annotations, pass `transfer="annotations"`:

```python
artifact = db.Artifact.get(key="example_datasets/mini_immuno/dataset1.h5ad")
artifact.save(transfer="annotations")
```

The artifact now has all feature & label annotations:

```python
artifact.describe()
```

The sync is zero-copy, which means that the data itself remained in the original storage location:

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

As expected, upon re-syncing an object, `lamindb` identifies that the object already exists in the target database and simply maps it:

```python
artifact = db.Artifact.get(key="example_datasets/mini_immuno/dataset1.h5ad")
artifact.save()
```

```{dropdown} How do I know if an object is in the default database or elsewhere?

Every `SQLRecord` object has an attribute `._state.db` which can take the following values:

- `None`: the object has not yet been saved to any database
- `"default"`: the object is saved on the default database instance
- `"account/name"`: the object is saved on a non-default database instance referenced by `account/name` (e.g., `laminlabs/lamindata`)

```

```python
# test the last 3 cells here
assert artifact.transform.description == "Transfer from `laminlabs/lamindata`"
assert artifact.transform.key == "__lamindb_transfer__/4XIuR0tvaiXM"
assert artifact.transform.uid == "4XIuR0tvaiXM0000"
assert artifact.run.initiated_by_run.transform.description.startswith("Sync data")
```
