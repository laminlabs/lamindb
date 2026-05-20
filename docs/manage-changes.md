# Manage changes

You can manage updates to objects through versioning, and safely remove obsolete data by moving it to the `archive` or `trash` branches. For more complex workflows, you can use contribution branches similar to branches in git.

## Versioning

## Archive & trash

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

If you delete an artifact, it gets moved onto the `trash` branch.

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

## Contribution branches

To create a contribution branch and switch to it, run:

```bash
lamin switch -c my_branch
```

This configures a default branch in your environment that takes effect in shell, Python, and R sessions. All objects that you create will then be created on that branch. Alternatively, you can directly configure a branch via the API:

::::{tab-set}

:::{tab-item} Via settings

```python
ln.setup.settings.branch = "my_branch"  # globally switch the default branch
```

:::

:::{tab-item} Via `ln.track()`

```python
ln.track(branch="my_branch")  # default branch for all objects created in a run to my_branch
```

:::

:::{tab-item} Via object constructor

```python
ln.Artifact(..., branch="my_branch")  # add an artifact on my_branch
ln.ULabel(..., branch="my_branch")  # add a ULabel on my_branch
```

:::
::::

To merge a contribution branch into `main`, run:

```bash
lamin switch main  # switch to the main branch
lamin merge my_branch  # merge contribution branch into main
```

To see the current branch along with other information, run:

```bash
lamin info
```

To annotate the current branch with a `README.md`, run:

```bash
lamin annotate branch --readme README.md
```

To comment on the current branch, run:

```bash
lamin annotate branch --comment "I think we should revisit this, tomorrow, WDYT?"
```

To describe the current branch (optionally include comments), run:

```bash
lamin describe branch --include comments
```

To trace on which branch a `SQLRecord` object was created, run:

```python
sqlrecord.created_on.describe()
```

To open a Change Request for a branch, run:

::::{tab-set}

:::{tab-item} CLI

```bash
lamin update branch --status draft  # for current branch
lamin update branch --name my_branch --status review  # for any branch
```

:::

:::{tab-item} Python

```python
branch = ln.Branch.get(name="my_branch")
branch.status = "draft"
branch.save()

branch.status = "review"
branch.save()
```

:::
::::

Just like Pull Requests on GitHub, branches are never deleted
so that the provenance of a change stays traceable.

For the API reference, see {class}`~lamindb.Branch`.
