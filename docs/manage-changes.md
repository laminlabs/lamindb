# Manage changes

To manage changes, you can use versioning, branching, an `archive` and the `trash`.

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

::::{tab-set}
:::{tab-item} UI

<img width="800" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/pQWzkbPHq1Sk2zs50001.png"/>
:::

:::{tab-item} CLI

```bash
lamin list branch
```

:::

:::{tab-item} Python

```python
ln.Branch.to_dataframe()
```

:::
::::

### Archive & trash

If you delete an object, it gets moved into the `trash`. There, it's hidden from queries and searches and scheduled for deletion.

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

Objects in the archive are hidden from queries and searches, like objects in the trash, but they are not scheduled for deletion. You can query for them by adding `branch__name="archive"` to the filter.

### Contribution workflow

#### Create a branch

To create a contribution branch and switch to it, run:

```bash
lamin switch -c my_branch
```

This configures a default branch in your environment that takes effect in shell, Python, and R sessions. All objects you create are then created on that branch.

:::::{dropdown} Alternatively, you can configure a branch via the API:

::::{tab-set}

:::{tab-item} Via settings

```python
ln.setup.settings.branch = "my_branch"  # equivalent to lamin switch my_branch via CLI
```

:::

:::{tab-item} Via `ln.track()`

```python
ln.track(branch="my_branch")  # default branch for all objects created in a run on my_branch
```

:::

:::{tab-item} Via constructor

```python
ln.Artifact(..., branch="my_branch")  # add an artifact on my_branch
ln.ULabel(..., branch="my_branch")  # add a ULabel on my_branch
```

:::
::::
:::::

#### Open a Change Request

::::{tab-set}

:::{tab-item} UI

Open the branch in the Changes page, and use the "Make change request" button to set it to "draft."

<img width="800" alt="branch make change request" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/HH5rRZrjrRw2cWa70000.png"/>

Once you think the branch is ready, use the "Mark ready for review" button to submit it for review.

<img width="800" alt="branch mark ready" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/WrZILsMtFsYfcc8V0000.png"/>

:::

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

#### Merge changes

::::{tab-set}

:::{tab-item} UI

To merge your contribution branch into the `main` branch, set the "Target Branch" dropdown to `main` and use the "Merge" button.

<img width="800" alt="branch merge changes" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/HOykpFPugonJyerM0000.png"/>

:::

:::{tab-item} CLI

```bash
lamin switch main  # switch to the main branch
lamin merge my_branch  # merge contribution branch into main
```

:::

:::{tab-item} Python

```python
ln.setup.merge("my_branch", target="main")
```

:::
::::

### Work with branches

To see the current branch along with other information, run:

```bash
lamin info
```

Add a branch README:
::::{tab-set}

:::{tab-item} UI
<img width="800" alt="branch readme" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/IIKfNn6m3vlAReXz0000.png"/>
:::

:::{tab-item} CLI
To annotate the current branch with a `README.md`, run:

```bash
lamin annotate branch --readme README.md
```

:::
::::

Comment on a branch:

::::{tab-set}

:::{tab-item} UI
<img width="800" alt="branch commenting" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/PNH3OlEeGLQ8sFpl0000.png"/>
:::

:::{tab-item} CLI
To comment on the current branch, run:

```bash
lamin annotate branch --comment "I think we should revisit this, tomorrow, WDYT?"
```

:::
::::

To describe the current branch (optionally include comments), run:

```bash
lamin describe branch --include comments
```

To trace on which branch a `SQLRecord` object was created, run:

```python
sqlrecord.created_on.describe()
```

Just like pull requests on GitHub, branches are never deleted so that the provenance of a change stays traceable.

For the API reference, see {class}`~lamindb.Branch`.
