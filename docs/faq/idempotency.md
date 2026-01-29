---
executable_via: python
---

# Will data get duplicated upon re-running code?

LaminDB's operations are idempotent in the sense defined here, which allows you to re-run code without duplicating data.

:::{admonition} SQLRecords with `name` field

When you instantiate {class}`~lamindb.models.SQLRecord` with a name, in case a name has an _exact match_ in a registry, the constructor returns it instead of creating a new record. In case records with _similar names_ exist, you'll see them in a table: you can then decide whether you want to save the new record or pick an existing record.

If you set {attr}`~lamindb.core.subsettings.CreationSettings.search_names` to `False`, you bypass these checks.

:::

:::{admonition} Artifacts & collections

If you instantiate {class}`~lamindb.Artifact` from data that already exists as an artifact, the `Artifact()` constructor returns the existing artifact based on a hash lookup.

:::

```python
# pip install lamindb
!lamin init --storage ./test-idempotency
```

```python
import lamindb as ln

ln.track("ANW20Fr4eZgM0000")
```

## SQLRecords with name field

```python
assert ln.settings.creation.search_names
```

Let us add a first record to the {class}`~lamindb.Record` registry:

```python
label = ln.Record(name="My label 1").save()
```

If we create a new record, we'll automatically get search results that give clues on whether we are prone to duplicating an entry:

```python
label = ln.Record(name="My label 1a")
```

Let's save the `1a` label, we actually intend to create it.

```python
label.save()
```

In case we match an existing name directly, we'll get the existing object:

```python
label = ln.Record(name="My label 1")
```

If we save it again, it will not create a new entry in the registry:

```python
label.save()
```

Now, if we create a third record, we'll get two alternatives:

```python
label = ln.Record(name="My label 1b")
```

If we prefer to not perform a search, e.g. for performance reasons, we can switch it off.

```python
ln.settings.creation.search_names = False
label = ln.Record(name="My label 1c")
```

Switch it back on:

```python
ln.settings.creation.search_names = True
```

## Artifacts & collections

```python
filepath = ln.examples.datasets.file_fcs()
```

Create an `Artifact`:

```python
artifact = ln.Artifact(filepath, key="my_fcs_file.fcs").save()
```

```python
assert artifact.hash == "rCPvmZB19xs4zHZ7p_-Wrg"
assert artifact.run == ln.context.run
assert not artifact.recreating_runs.exists()
```

Create an `Artifact` from the same path:

```python
artifact2 = ln.Artifact(filepath, key="my_fcs_file.fcs")
```

It gives us the existing object:

```python
assert artifact.id == artifact2.id
assert artifact.run == artifact2.run
assert not artifact.recreating_runs.exists()
```

If you save it again, nothing will happen (the operation is idempotent):

```python
artifact2.save()
```

In the hidden cell below, you'll see how this interplays with data lineage.

```python
ln.track(new_run=True)
artifact3 = ln.Artifact(filepath, key="my_fcs_file.fcs")
assert artifact3.id == artifact2.id
assert artifact3.run == artifact2.run != ln.context.run  # run is not updated
assert artifact2.recreating_runs.first() == ln.context.run
```

```python
!rm -rf ./test-idempotency
!lamin delete --force test-idempotency
```
