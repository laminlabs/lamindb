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

If you already have a set of artifacts and you'd like to stream their content, see {doc}`arrays`.

```python
# initialize a test database to run examples
!lamin init --storage ./test-registries --modules bionty
```

Let's start by creating a few exemplary datasets:

```python
import lamindb as ln

ln.Artifact(ln.examples.datasets.file_fastq(), key="raw/my_fastq.fastq.gz").save()
ln.Artifact(ln.examples.datasets.file_jpg_paradisi05(), key="my_image.jpg").save()
ln.Artifact.from_dataframe(ln.examples.datasets.df_iris(), key="iris.parquet").save()
ln.examples.datasets.mini_immuno.save_mini_immuno_datasets()
```

## Get an overview

The easiest way to get an overview over all artifacts is by typing {meth}`~lamindb.Artifact.to_dataframe`, which returns the most recently created artifacts in the {class}`~lamindb.Artifact` registry.

```python
ln.Artifact.to_dataframe()
```

You can include features.

```python
ln.Artifact.to_dataframe(include="features")
```

You can include fields from other registries.

```python
ln.Artifact.to_dataframe(
    include=[
        "created_by__name",
        "records__name",
        "cell_types__name",
        "schemas__itype",
    ]
)
```

You can also get an overview of the entire database.

```python
ln.view()
```

## Auto-complete objects

For registries with less than 100k objects, auto-completing a `Lookup` object is the most convenient way of finding a record.

```python
records = ln.Record.lookup()
```

With auto-complete, we find a record:

```python
experiment_1 = records.experiment_1
experiment_1
```

This works for any {class}`~lamindb.models.BaseSQLRecord` class, e.g., also for plugin `bionty`.

```python
import bionty as bt

cell_types = bt.CellType.lookup()
```

## Get

{meth}`~lamindb.models.BaseSQLRecord.get` errors if none or more than one matching objects are found.

```python
ln.Record.get(experiment_1.uid)  # by uid
ln.Record.get(name="Experiment 1")  # by field
```

## Search

You can search every registry via {meth}`~lamindb.models.BaseSQLRecord.search`. For example, the `Artifact` registry.

```python
ln.Artifact.search("iris").to_dataframe()
```

Here is more background on search and examples for searching the entire cell type ontology: {doc}`/faq/search`

## Queries

### By fields

Use {meth}`~lamindb.models.BaseSQLRecord.filter` to query all artifacts by the `suffix` field:

```python
qs = ln.Artifact.filter(suffix=".h5ad")
qs
```

This returns a {class}`~lamindb.models.QuerySet`, which lazily references the set of {class}`~lamindb.models.BaseSQLRecord` objects that matches the filter statement. You can iteratively filter a queryset:

```python
qs = qs.filter(records__name="Experiment 1")
```

To access the results encoded in a queryset, call:

- {meth}`~lamindb.models.BasicQuerySet.to_dataframe`: A pandas `DataFrame` with each record in a row.
- {meth}`~lamindb.models.BasicQuerySet.one`: Exactly one record. Will raise an error if there is none. Is equivalent to the `.get()` method shown above.
- {meth}`~lamindb.models.BasicQuerySet.one_or_none`: Either one record or `None` if there is no query result.

Alternatively,

- use the `QuerySet` as an iterator
- get individual objects via `qs[0]`, `qs[1]`

For example:

```python
qs.to_dataframe()
```

The `SQLRecord` classes in LaminDB are Django Models and any [Django query](https://docs.djangoproject.com/en/stable/topics/db/queries/) works.

Django has a double-under-score syntax to filter based on related tables.
This syntax enables you to traverse several layers of relations and comparators.

```python
ln.Artifact.filter(created_by__handle__startswith="testuse").to_dataframe()
```

The filter selects all artifacts based on the users who ran the generating notebook. Under the hood, in the SQL database, it's joining the artifact table with the user table.

Another example would be querying all datasets that measure a particular feature. For instance, which datasets measures `"CD8A"`. Here is how:

```python
cd8a = bt.Gene.get(symbol="CD8A")
# query for all feature sets that contain CD8A
schemas_with_cd8a = ln.Schema.filter(genes=cd8a)
# get all artifacts
ln.Artifact.filter(schemas__in=schemas_with_cd8a).to_dataframe()
```

Instead of splitting this across three queries, the double-underscore syntax allows you to define a path for one query.

```python
ln.Artifact.filter(schemas__genes__symbol="CD8A").to_dataframe()
```

### By features

The `Artifact`, `Record`, and `Run` registries can be queried by features, via an implicit lookup in the {class}`~lamindb.Feature` registry:

<!-- #region -->
<!-- cannot run tabbed code on CI, see test_artifact_filter_by_multiple_features for tests -->

::::{tab-set}

:::{tab-item} Via strings / kwargs

```python
ln.Artifact.filter(
    perturbation="DMSO",
    temperature__gt=26,
).to_dataframe(include="features")
```

:::

:::{tab-item} Via objects / expressions

```python
perturbation = ln.Feature.get(name="perturbation")  # can optionally pass a feature type to disambiguate
temperature = ln.Feature.get(name="temperature")
ln.Artifact.filter(    # note this is now an expression using the == syntax
    perturbation == "DMSO",
    temperature > 21,
).to_dataframe(include="features")
```

:::

::::

<!-- #endregion -->

Just like for fields holding dictionary values, you can query for dictionary keys in features whose `dtype` is `dict`:

```python
ln.Artifact.filter(study_metadata__detail1="123").to_dataframe(include="features")
```

You can query for whether a dataset is annotated annotated by a feature:

```python
ln.Artifact.filter(perturbation__isnull=False).to_dataframe(include="features")
```

## Filter operators

You can qualify the type of comparison in a query by using a comparator.

Below follows a list of the most import, but Django supports about [two dozen field comparators](https://docs.djangoproject.com/en/stable/ref/models/querysets/#field-lookups) `field__comparator=value`.

### and

```python
ln.Artifact.filter(suffix=".h5ad", records=experiment_1).to_dataframe()
```

### less than/ greater than

Or subset to artifacts greater than 10kB. Here, we can't use keyword arguments, but need an explicit where statement.

```python
ln.Artifact.filter(records=experiment_1, size__gt=1e4).to_dataframe()
```

### in

```python
ln.Artifact.filter(suffix__in=[".jpg", ".fastq.gz"]).to_dataframe()
```

### order by

```python
ln.Artifact.filter().order_by("created_at").to_dataframe()
```

```python
# reverse ordering
ln.Artifact.filter().order_by("-created_at").to_dataframe()
```

```python
ln.Artifact.filter().order_by("key").to_dataframe()
```

```python
# reverse ordering
ln.Artifact.filter().order_by("-key").to_dataframe()
```

### contains

```python
ln.Transform.filter(description__contains="search").to_dataframe().head(5)
```

And case-insensitive:

```python
ln.Transform.filter(description__icontains="Search").to_dataframe().head(5)
```

### startswith

```python
ln.Transform.filter(description__startswith="Query").to_dataframe()
```

### or

```python
ln.Artifact.filter(ln.Q(suffix=".jpg") | ln.Q(suffix=".fastq.gz")).to_dataframe()
```

### negate/ unequal

```python
ln.Artifact.filter(~ln.Q(suffix=".jpg")).to_dataframe()
```

### JSON

Here is an example for querying by parameters: {ref}`track-run-parameters`.
