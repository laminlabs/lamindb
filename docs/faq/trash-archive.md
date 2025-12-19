# How do I trash or archive objects?

Any object in LaminDB has the following 3 levels of visibility through 3 default branches:

- `main`: visible
- `archive`: excluded from query & search
- `trash`: excluded from query & search, scheduled for deletion

Let's look at an example for an `Artifact` object while noting that the same applies to any other `SQLRecord`.

```python
import lamindb as ln
import pandas as pd

df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
artifact = ln.Artifact.from_dataframe(df, key="dataset.parquet").save()
```

An artifact is by default created on the `main` branch.

```python
assert artifact.branch.name == "main"
ln.Artifact.filter(key="dataset.parquet").to_dataframe()
# the artifact shows up
```

If you delete an artifact, it gets moved into the `trash` branch.

```python
artifact.delete()
assert artifact.branch.name == "trash"
```

Artifacts in trash won't show up in queries with default arguments:

```python
ln.Artifact.filter(key="dataset.parquet").to_dataframe()
# the artifact does not show up
```

You can query for them by adding the `trash` branch to the filter.

```python
ln.Artifact.filter(key="dataset.parquet", branch__name="trash").to_dataframe()
# the artifact shows up
```

You can restore an artifact from trash:

```python
artifact.restore()
ln.Artifact.filter(key="dataset.parquet").to_dataframe()
# the artifact shows up
```
