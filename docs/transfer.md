---
executable_via: python
---

# Transfer data

This guide shows how to transfer data from a source database into the currently connected database.

```python
# pip install lamindb
!lamin init --storage ./test-transfer --modules bionty
```

```python
import lamindb as ln

ln.track()
```

Query all artifacts in the `laminlabs/lamindata` instance and filter them to their latest versions.

```python
# query all latest artifact versions
artifacts = ln.Artifact.connect("laminlabs/lamindata").filter(is_latest=True)

# convert the QuerySet to a DataFrame and show the latest 5 versions
artifacts.to_dataframe().head()
```

You can now further subset or search the {class}`~lamindb.models.QuerySet`. Here we query by whether the description contains "tabula sapiens".

```python
artifact = artifacts.filter(description__contains="Tabula Sapiens").first()
artifact.describe()
```

By saving the artifact record that's currently attached to the source database instance, you transfer it to the default database instance.

```python
artifact.save()
```

```{dropdown} How do I know if a record is saved in the default database instance or not?

Every record has an attribute `._state.db` which can take the following values:

- `None`: the record has not yet been saved to any database
- `"default"`: the record is saved on the default database instance
- `"account/name"`: the record is saved on a non-default database instance referenced by `account/name` (e.g., `laminlabs/lamindata`)

```

The artifact record has been transferred to the current database without feature & label annotations, but with updated data lineage.

```python
artifact.describe()
```

You see that the data itself remained in the original storage location, which has been added to the current instance's storage location as a read-only location (indicated by the fact that the `instance_uid` doesn't match the current instance).

```python
ln.Storage.to_dataframe()
```

See the state of the database.

```python
ln.view()
```

View lineage:

```python
artifact.view_lineage()
```

The transferred dataset is linked to a special type of transform that stores the slug and uid of the source instance:

```python
artifact.transform.description
```

The transform key has the form `f"__lamindb_transfer__/{source_instance.uid}"`:

```python
artifact.transform.key
```

The current notebook run is linked as the initiated_by_run of the "transfer run":

```python
artifact.run.initiated_by_run.transform
```

Upon re-transferring a record, it will identify that the record already exists in the target database and simply map the record.

```python
artifact = artifacts.filter(description__contains="Tabula Sapiens").first()
artifact.save()
```

If you also want to transfer annotations of the artifact, you can pass `transfer="annotations"` to `save()`. Just note that this might populate your target database with metadata that doesn't match the conventions you want to enforce.

```python
artifact = artifacts.filter(description__contains="Tabula Sapiens").first()
artifact.save(transfer="annotations")
```

The artifact is now annotated.

```python
artifact.describe()
```

```python
# test the last 3 cells here
assert artifact.transform.description == "Transfer from `laminlabs/lamindata`"
assert artifact.transform.key == "__lamindb_transfer__/4XIuR0tvaiXM"
assert artifact.transform.uid == "4XIuR0tvaiXM0000"
assert artifact.run.initiated_by_run.transform.description == "Transfer data"
```
