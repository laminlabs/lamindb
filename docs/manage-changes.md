# Manage changes

Managing changes in LaminDB is largely analogous to managing code changes via branching in git and Pull Requests in GitHub.

To create a contribution {class}`~lamindb.Branch` and switch to it, run:

```bash
lamin switch -c my_branch
```

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
