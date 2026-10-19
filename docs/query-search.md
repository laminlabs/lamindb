# Query, search & load data

If you just want to use `lamindb` like a file system, then you can do that by querying artifacts through their file keys.

```python
import lamindb as ln

# one file, latest version
artifact = ln.Artifact.get(key="project1/dataset1.csv")

# a previous version of the same file
artifact_v1 = ln.Artifact.get(key="project1/dataset1.csv", version="1")

# every file in a folder
artifacts = ln.Artifact.filter(key__startswith="project1/")

# load into memory, or sync to a local path
df = artifact.load()
path = artifact.cache()
```

If you're looking for a data catalog, a data lake, or a lakehouse, see the following guides.

```{toctree}
:maxdepth: 1
:hidden:

registries
tables
arrays
```
