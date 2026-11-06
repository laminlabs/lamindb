# Query, search & load data

If you just want to use `lamindb` like a file system, then you can do that by querying artifacts through their file keys.

```python
import lamindb as ln

# one file, latest version
artifact = ln.Artifact.get(key="project1/dataset1.csv")

# a previous version of the same file
artifact_v1 = ln.Artifact.get(key="project1/dataset1.csv", version="1")

# sync file to a local cache
path = artifact.cache()

# a folder stored as one artifact
folder = ln.Artifact.get(key="project1/my_folder")

# sync the folder into a local cache
folder_path = folder.cache()

# all artifacts whose file key starts with "project1/"
artifacts = ln.Artifact.filter(key__startswith="project1/")
```

## Aggregate & annotate

`QuerySet` supports Django's `annotate()` for computing aggregates over related records and filtering on them:

```python
from django.db.models import Count

# artifacts by the number of schemas they validate against
ln.Artifact.filter().annotate(n_schemas=Count("schemas")).filter(n_schemas=2)
```

You can traverse any relationship defined on your registries. For example, if your instance links artifacts to a custom `GeneticTreatment` registry, count the target genes per artifact:

```python
from django.db.models import Count

ln.Artifact.filter().annotate(
    n_genes=Count("genetic_treatments__targets__genes")
).filter(n_genes=4)
```

If you're looking for a data catalog, a data lake, or a lakehouse, see the following guides:

```{toctree}
:maxdepth: 1

registries
tables
arrays
```
