---
execute_via: python
---

# Query & search

```{toctree}
:maxdepth: 1
:hidden:

parquet
arrays
```

This guide walks through querying & searching registries.
To understand the links between objects that enable relational queries, see {doc}`organize`.
To stream datasets from storage, see {doc}`arrays`.

## Basics

Create a database object:

```python
import lamindb as ln

db = ln.DB("laminlabs/lamindata")
```

### Get an overview

Return the most recent entries of the {class}`~lamindb.Artifact` registry as a `DataFrame`:

```python
db.Artifact.to_dataframe()
```

Include fields from other registries with the `__` syntax:

```python
db.Artifact.to_dataframe(include=["created_by__handle"])
```

Get an overview of all most recent objects in the database:

```python
db.view()
```

### Get one object

{meth}`~lamindb.models.BaseSQLRecord.get` errors if it doesn't find exactly one matching object:

```python
db.Record.get("HhKWgzjVW6cfnVep")  # by uid
db.Record.get(name="EXP-RNA-001")  # by field, here the name field
```

### Search

You can search every registry via {meth}`~lamindb.models.BaseSQLRecord.search`. For example, the `Artifact` registry:

```python
db.Artifact.search("iris").to_dataframe()
```

Here is more background on search and examples for searching the cell type registry: {doc}`/faq/search`

Auto-completing can also be a good way of finding an object:

```python
records = db.Record.lookup()
experiment_1 = records.exp_rna_001  # assuming a record named 'EXP-RNA-001' exists
experiment_1
```

This works for any {class}`~lamindb.models.BaseSQLRecord` class, e.g., also for plugin `bionty`:

```python
cell_types = db.bionty.CellType.lookup()
```

## Queries

### By fields

Use {meth}`~lamindb.models.BaseSQLRecord.filter` to query artifacts by fields, e.g., by {attr}`~lamindb.Artifact.suffix`:

```python
qs = db.Artifact.filter(suffix=".csv")
qs
```

This returns a {class}`~lamindb.models.QuerySet` with the objects that match the filter. You can filter a queryset:

```python
qs = qs.filter(records=experiment_1)  # filter to artifacts annotated with a record
qs.to_dataframe()
```

To access the results in a queryset:

- {meth}`~lamindb.models.BasicQuerySet.to_dataframe`: A pandas `DataFrame` with each object in a row.
- {meth}`~lamindb.models.BasicQuerySet.one`: Exactly one object. Will raise an error if there is none. Is equivalent to the `.get()` method shown above.
- {meth}`~lamindb.models.BasicQuerySet.one_or_none`: Either one object or `None` if there is no query result.

You can also:

- use the `QuerySet` as an iterator
- get individual objects via `qs[0]`, `qs[1]`

The `SQLRecord` classes (registries) in LaminDB are Django Models and any [Django query](https://docs.djangoproject.com/en/stable/topics/db/queries/) works.
Django has a double-under-score syntax to filter based on related tables.
This syntax enables you to traverse several layers of relations and comparators.

For example, the following filter selects artifacts based on the users who ran the generating notebook. Under the hood, in the SQL database, it's joining the artifact table with the user table.

```python
db.Artifact.filter(created_by__handle__startswith="testuse").to_dataframe()
```

Another example would be querying datasets that measure a particular feature. For instance, which datasets measure expression of `CD8A`:

```python
cd8a = db.bionty.Gene.get(symbol="CD8A")
# query for all schemas that contain CD8A
schemas_with_cd8a = db.Schema.filter(genes=cd8a)
# get all artifacts
db.Artifact.filter(schemas__in=schemas_with_cd8a).to_dataframe()
```

Instead of splitting this across three queries, the double-underscore syntax allows you to define a path for one query:

```python
db.Artifact.filter(schemas__genes__symbol="CD8A").to_dataframe()
```

### By features

The {class}`~lamindb.Feature` registry indexes variables across datasets to enable querying by dimensions. Registry fields couldn't scale to millions of dimensions as found in biology:

<img width="800px" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/VFFgFdAlJnssyOdk0001.svg">

The {class}`~lamindb.Artifact`, {class}`~lamindb.Record`, and {class}`~lamindb.Run` registries can be queried by feature expressions:

```python
feature_type = db.Feature.get(name="mini_immuno")
perturbation = db.Feature.get(name="perturbation", type=feature_type)
temperature = db.Feature.get(name="temperature")
db.Artifact.filter(
    perturbation == "DMSO",
    temperature > 21,
).to_dataframe(include="features")
```

You can query for whether a dataset is annotated by a feature:

```python
db.Artifact.filter(perturbation.is_null(False)).to_dataframe(include="features")
```

## Cheat sheet: comparators

You can qualify the type of comparison in a query by using a comparator.
Below is a list of the most important ones, but Django supports about [two dozen field comparators](https://docs.djangoproject.com/en/stable/ref/models/querysets/#field-lookups) `field__comparator=value`.

### and

```python
db.Artifact.filter(suffix=".h5ad", records=experiment_1).to_dataframe()
```

### less than/ greater than

```python
# artifacts greater than 10kB
db.Artifact.filter(size__gt=1e4).to_dataframe()
```

### in

```python
db.Artifact.filter(suffix__in=[".jpg", ".fastq.gz"]).to_dataframe()
```

### order by

```python
db.Artifact.filter().order_by("created_at").to_dataframe()
```

```python
# reverse ordering
db.Artifact.filter().order_by("-created_at").to_dataframe()
```

```python
db.Artifact.filter().order_by("key").to_dataframe()
```

```python
# reverse ordering
db.Artifact.filter().order_by("-key").to_dataframe()
```

### contains

```python
db.Transform.filter(description__contains="search").to_dataframe().head(5)
```

And case-insensitive:

```python
db.Transform.filter(description__icontains="Search").to_dataframe().head(5)
```

### startswith

```python
db.Transform.filter(description__startswith="Query").to_dataframe()
```

### or

LaminDB uses Django's `Q` objects to encapsulate filters for complex boolean logic.

```python
db.Artifact.filter(ln.Q(suffix=".jpg") | ln.Q(suffix=".fastq.gz")).to_dataframe()
```

### negate/ unequal

```python
db.Artifact.filter(~ln.Q(suffix=".jpg")).to_dataframe()
```

### JSON

```python
db.Run.filter(params__learning_rate__gt=0.01).to_dataframe()
```

Here is an example for querying by parameters: {ref}`track-run-parameters`.
