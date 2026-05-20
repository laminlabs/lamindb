# Manage changes

You can manage updates to objects through versioning, and safely remove obsolete data by moving it to the `archive` or `trash` branches. For more complex workflows, you can use contribution branches similar to branches in git.

## Versioning

You can make a new version of an artifact, transform, or collection by passing an existing `key`. For example, for an artifact:

```python
import lamindb as ln
from pathlib import Path

Path("my_file.txt").write_text("v1")
artifact = ln.Artifact("my_file.txt", key="my_file.txt").save()

Path("my_file.txt").write_text("v2")
artifact_v2 = ln.Artifact("my_file.txt", key="my_file.txt").save()

artifact_v2.versions.to_dataframe()  # see all versions of this artifact
```

This works because {class}`~lamindb.Artifact`, {class}`~lamindb.Transform`, and {class}`~lamindb.Collection` inherit from {class}`~lamindb.models.IsVersioned`.

## Branching

All primary objects like artifacts, records, transforms, etc. (any that inherit from {class}`~lamindb.models.SQLRecord`) have a {attr}`~lamindb.models.SQLRecord.branch` field that determines their life cycle.

There are three built-in branches: `main`, `trash`, and `archive`. By default, objects are created on the `main` branch and visible in queries and searches.

### Archive & trash

If you delete an object, it gets moved into the `trash`. There, it's hidden from queries & search and scheduled for deletion.

```python
artifact.delete()
assert artifact.branch.name == "trash"

# the artifact does not show up in a default query
ln.Artifact.filter(key="my_file.txt")

# you can still query for it by adding the trash branch to the filter
ln.Artifact.filter(key="my_file.txt", branch__name="trash")

# you can restore it from trash
artifact.restore()
```

To move an object into the archive, run:

```python
artifact.branch_id = 0
artifact.save()
```

Objects in the archive are hidden from queries & search like objects in the trash, but they are not scheduled for deletion. You can query for them by adding `branch__name="archive"` to the filter.

### Contribution branches

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
