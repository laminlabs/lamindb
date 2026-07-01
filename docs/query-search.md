---
execute_via: python
---

# Query & search

```{toctree}
:maxdepth: 1
:hidden:

arrays
```

This guide walks through different ways of querying & searching registries.
To understand the underlying cross-linking of objects in the SQL database, see {doc}`organize`.
To stream datasets from storage, see {doc}`arrays`.

## Basics

Create a database object:

```python
import lamindb as ln

db = ln.DB("laminlabs/lamindata")
```

### Get an overview

The easiest way to get an overview over all artifacts is by typing {meth}`~lamindb.Artifact.to_dataframe`, which returns the most recent entries of the {class}`~lamindb.Artifact` registry:

```python
db.Artifact.to_dataframe()
```

You can include fields from other registries with the `__` syntax:

```python
db.Artifact.to_dataframe(include=["created_by__handle"])
```

Get an overview of the most recent objects in the database:

```python
db.view()
```

### Auto-complete

For registries with less than 100k objects, auto-completing a `Lookup` object is a good way of finding an object.

```python
records = db.Record.lookup()
```

With auto-complete, we find an object:

```python
experiment_1 = records.experiment_1
experiment_1
```

This works for any {class}`~lamindb.models.BaseSQLRecord` class, e.g., also for plugin `bionty`.

```python
cell_types = db.bionty.CellType.lookup()
```

### Get one object

{meth}`~lamindb.models.BaseSQLRecord.get` errors if none or more than one matching objects are found.

```python
db.Record.get(experiment_1.uid)  # by uid
db.Record.get(name="Experiment 1")  # by field
```

### Search

You can search every registry via {meth}`~lamindb.models.BaseSQLRecord.search`. For example, the `Artifact` registry.

```python
db.Artifact.search("iris").to_dataframe()
```

Here is more background on search and examples for searching the entire cell type ontology: {doc}`/faq/search`

## Queries

### By fields

Use {meth}`~lamindb.models.BaseSQLRecord.filter` to query artifacts by any field, e.g., the {attr}`~lamindb.Artifact.suffix` field:

```python
qs = db.Artifact.filter(suffix=".h5ad")
qs
```

This returns a {class}`~lamindb.models.QuerySet`, which lazily references the set of {class}`~lamindb.models.BaseSQLRecord` objects that matches the filter. You can filter a queryset:

```python
qs = qs.filter(records__name="Experiment 1")
qs.to_dataframe()
```

To access the results encoded in a queryset, you can call one of:

- {meth}`~lamindb.models.BasicQuerySet.to_dataframe`: A pandas `DataFrame` with each record in a row.
- {meth}`~lamindb.models.BasicQuerySet.one`: Exactly one record. Will raise an error if there is none. Is equivalent to the `.get()` method shown above.
- {meth}`~lamindb.models.BasicQuerySet.one_or_none`: Either one record or `None` if there is no query result.

Alternatively,

- use the `QuerySet` as an iterator
- get individual objects via `qs[0]`, `qs[1]`

The `SQLRecord` classes in LaminDB are Django Models and any [Django query](https://docs.djangoproject.com/en/stable/topics/db/queries/) works.
Django has a double-under-score syntax to filter based on related tables.
This syntax enables you to traverse several layers of relations and comparators.

For example, the following filter selects artifacts based on the users who ran the generating notebook. Under the hood, in the SQL database, it's joining the artifact table with the user table.

```python
db.Artifact.filter(created_by__handle__startswith="testuse").to_dataframe()
```

Another example would be querying datasets that measure a particular feature. For instance, which datasets measures expression of `CD8A`:

```python
cd8a = db.bionty.Gene.get(symbol="CD8A")
# query for all feature sets that contain CD8A
schemas_with_cd8a = db.Schema.filter(genes=cd8a)
# get all artifacts
db.Artifact.filter(schemas__in=schemas_with_cd8a).to_dataframe()
```

Instead of splitting this across three queries, the double-underscore syntax allows you to define a path for one query:

```python
db.Artifact.filter(schemas__genes__symbol="CD8A").to_dataframe()
```

### By features

The {class}`~lamindb.Feature` registry indexes variables across datasets to enable querying by dimensions. Registry fields couldn't scale to millions of dimensions.

<img width="800px" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/VFFgFdAlJnssyOdk0001.svg">

The {class}`~lamindb.Artifact`, {class}`~lamindb.Record`, and {class}`~lamindb.Run` registries can be queried by features:

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
Below follows a list of the most important, but Django supports about [two dozen field comparators](https://docs.djangoproject.com/en/stable/ref/models/querysets/#field-lookups) `field__comparator=value`.

### and

```python
db.Artifact.filter(suffix=".h5ad", records=experiment_1).to_dataframe()
```

### less than/ greater than

```python
# artifacts greater than 10kB
db.Artifact.filter(records=experiment_1, size__gt=1e4).to_dataframe()
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

```python
db.Artifact.filter(ln.Q(suffix=".jpg") | ln.Q(suffix=".fastq.gz")).to_dataframe()
```

### negate/ unequal

```python
db.Artifact.filter(~ln.Q(suffix=".jpg")).to_dataframe()
```

### JSON

Here is an example for querying by parameters: {ref}`track-run-parameters`.
