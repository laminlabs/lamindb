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

::::{tab-set}
:::{tab-item} UI 
LaminHub has a Changes tab, which allows the user to have an overview of all the branches that have been created in the instance. Clicking on an individual branch page allows the user to view the entities like artifacts, transforms,etc. corresponding to the branch. 

<img width="400" alt="branch list page" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/pQWzkbPHq1Sk2zs50000.png"/>
:::

:::{tab-item} CLI
List all the branches in an instance by running

```bash
lamin list branch
```

:::
::::

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

To open a Change Request for a branch:

::::{tab-set}

:::{tab-item} UI

Open the branch in the Changes page, and use the 'Make change request' button to make it a 'draft' branch.
<img width="400" alt="branch make change request" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/HH5rRZrjrRw2cWa70000.png"/>

Once you think the branch is ready to be reviewed, use the 'Mark ready for review' button to submit it for review.
<img width="400" alt="branch mark ready" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/WrZILsMtFsYfcc8V0000.png"/>

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

Merging a contribution branch into `main`:

::::{tab-set}

:::{tab-item} UI

Use the 'Target Branch' dropdown to select 'main' as the target branch for the merge and use the 'Merge' button. 
<img width="400" alt="branch merge changes" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/HOykpFPugonJyerM0000.png"/>

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

To see the current branch along with other information, run:

```bash
lamin info
```

Adding a branch readme
::::{tab-set}

:::{tab-item} UI
<img width="400" alt="branch readme" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/IIKfNn6m3vlAReXz0000.png"/>
:::

:::{tab-item} CLI
To annotate the current branch with a `README.md`, run:
```bash
lamin annotate branch --readme README.md
```
:::
::::

Commenting on a branch

::::{tab-set}

:::{tab-item} UI
<img width="400" alt="branch commenting" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/PNH3OlEeGLQ8sFpl0000.png"/>
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


Just like Pull Requests on GitHub, branches are never deleted
so that the provenance of a change stays traceable.

For the API reference, see {class}`~lamindb.Branch`.
