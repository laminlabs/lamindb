# Organize datasets

```{toctree}
:maxdepth: 1
:hidden:

curate
```

This guide walks through organizing datasets using files & folders, database relationships, and versioned collections.

## Via files & folders

You can use LaminDB like a file system. Similar to AWS S3, you organize artifacts into virtual folders using `/`-separated keys. To ingest a single file into a `project1/` folder, you'd call:

```python
artifact1 = ln.Artifact("./dataset.csv", key="project1/dataset1.csv").save()
```

For convenience, if you want to create an artifact for every file in a directory, use {meth}`~lamindb.Artifact.from_dir`:

```python
artifacts = ln.Artifact.from_dir("./project1/").save()
```

You can then query for all artifacts in the `"./project1/"` folder via:

```python
artifacts = ln.Artifact.filter(key__startswith="project1/")
```

Unlike a regular file system, every artifact is versioned and comes with rich metadata.

:::{dropdown} What if I do not care about the metadata and version of every file in a folder?

In some cases a folder _is_ the dataset and you don't need fine-grained information for every file.
In this scenario, save the entire directory as a single artifact:

```python
ln.Artifact("./folder_abc", key="folder_abc").save()
```

:::

## Via relationships in the database

### Annotating with projects

What if an artifact is relevant to multiple projects?
A dataset that's in the `project1/` folder cannot **also** reside in a `project2/` folder.
You can solve this problem with the `artifact.projects` relationship that links the {class}`~lamindb.Project` to {class}`~lamindb.Artifact`:

<img width="400" alt="image" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/uVm5ptyqukPEKCix0000.png"/>

Here is how to annotate one artifact with two projects:

```python
project1 = ln.Project(name="Project 1").save()  # create project 1
project2 = ln.Project(name="Project 2").save()  # create project 2
artifact1.projects.add(project1, project2)      # annotate artifact1
```

This allows you to retrieve `artifact1` by querying any project it belongs to:

```python
artifacts_in_project1 = ln.Artifact.filter(projects=project1)
artifacts_in_project2 = ln.Artifact.filter(projects=project2)
```

Here, `artifact1` is part of both query results.

:::{dropdown} Three additional advantages to using related registries rather than folder structures.

1. Projects can be richly annotated (e.g., with start/end dates, parent projects, or member roles).
2. You no longer need to rely on fragile file paths. If a folder is renamed, path-based retrieval breaks, but a project query by `uid` will always work.[^protectproject]
3. You can run a constrained query or search against all projects in your database rather than trying to narrow a search to folder names.

:::

### Annotating with labels

You can annotate with other entity types, not just projects. LaminDB offers two main classes for this: {class}`~lamindb.Record` for metadata records and {class}`~lamindb.ULabel` for simple labels, which are both link to artifacts:

<img width="400" alt="image" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/qvhxt6UuoUO2Bd820000.png"/>

Here is how to annotate with a ulabel and with a sample record:

```python
ulabel1 = ln.ULabel(name="raw_data").save()  # create a ulabel
artifact1.ulabels.add(ulabel1)               # annotate artifact1

sample_type = ln.Record(                     # create a record type "Samples"
    name="Samples",
    is_type=True
).save()
record1 = ln.Record(                         # create a sample record
    name="My sample",
    features={"gc_content": 0.5}
).save()
artifact1.records.add(record1)               # annnotate artifact1
```

You can use records and ulabels alongside entity types in modules such as {mod}`bionty`:

```python
import bionty as bt

cell_type1 = bt.CellType.from_source(
    name="T cell"                            # create a cell type from a public ontology
).save()
artifact1.cell_types.add(cell_type1)         # annotate artifact1
```

### Annotating with features

To annotate with non-categorical data types or to disambiguate categorical annotations, use {class}`~lamindb.Feature` objects.

<img width="400" alt="image" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/eT6SEny5HpQQNgFl0000.png"/>

Here is how to define features and annotate an artifact with feature values:

```python
exp_type = ln.Record.get(name="Experiments")          # query the entity type `Experiments`
ln.Feature(name="gc_content", dtype=float).save()     # define a feature with dtype float
ln.Feature(name="experiment", dtype=exp_type).save()  # define a feature with dtype `Experiments`
artifact.features.set_values({
    "gc_content": 0.55,                               # validated to be a float
    "experiment": "Experiment 1",                     # validated to exist under the `Experiments` record type
})
```

When you work with structured data formats like `DataFrame` or `AnnData`, it often makes sense to validate the content of their features. After validation, the parsed feature values are automatically used for annotation. The easiest way is to use validation and auto-annotation is the built-in schema `"valid_features"`:

```python
# validate columns in the dataframe and map them on features
# auto-annotate with parsed metadata
ln.Artifact.from_dataframe(df, schema="valid_features").save()
```

Below is an example from the {doc}`docs:tutorial` illustrating how you get e.g. cell type, treatment, and assay annotations based on a dataframe's content. You can read more on this in {doc}`/curate`.

<img width="600px" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/6sofuDVvTANB0f480003.png">

### Annotating with data-lineage

When you call {func}`~lamindb.track` or decorate a function with {func}`~lamindb.flow`, you automatically annotate artifacts with {class}`~lamindb.Run` and {class}`~lamindb.Transform` objects.

<img width="400" alt="image" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/Z1iliqp5mInQQ2iY0000.png"/>

Here is how:

```{eval-rst}
.. literalinclude:: scripts/run_track_and_finish.py
   :language: python
```

Note that you can pass `project` to {func}`~lamindb.track` to auto-annotate all objects that are created in a run with a project label. Read more in {doc}`/track`.

### Overview of auto-generated annotations

The {class}`~lamindb.Artifact` registry has simple fields (such as `description`, `created_at`, `size`) and related fields (such as `projects`, `created_by`, `storage`). Many of these fields are automatically populated and you can use them to retrieve sets of artifacts.

<img width="800px" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/HMfWLa1rFkxcxQEN0000.svg">

All other registries link to {class}`~lamindb.Artifact` to provide context for finding, querying, validating, and managing artifacts.[^starsnowflake]

:::{dropdown} Can you give me some example queries?

Here are examples leveraging auto-populated fields.

```python
artifacts = ln.Artifact.filter(
    created_at__gt="2023-06-24",    # created after June 24th, 2023
    size__lt=1e9,                   # smaller than 1GB
    suffix=".parquet",              # with a .parquet suffix
    n_observations__gt=1000,        # with more than 1000 observations
    n_files__gt=1000,               # folder-like artifacts with more than 1000 files
    otype="DataFrame",              # that are DataFrames
    created_on__name="my-branch",   # created on a specific branch or environment
    created_by__handle="falexwolf", # created by user with handle falexwolf
    run=run,                        # created by a specific run
    transform__name="my-script.py", # created by a specific script/notebook
)
```

:::

## Versioned collections of artifacts

Sometimes, you need to both group artifacts by metadata and version the entire set. For this, use {class}`~lamindb.Collection`

<img width="160" alt="image" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/QR0KuktVEnVL08K90000.png"/>

Unlike during annotation, you have to pass an entire group of artifacts to a `Collection` constructor:

```python
collection = ln.Collection([artifact1, artifact2], key="my_data_release").save()
```

And unlike the folder-based or annotation-based sets of artifacts — which can change as artifacts are added or removed — a collection guarantees an exact, immutable set of artifacts.

Artifacts are versioned based on the hash of their content. Collections are versioned based on the top-level hash of their artifact hashes. If you use the {meth}`~lamindb.Collection.append` method, a new version of the collection is created, and the old version is left unchanged:

```python
collection_v2 = collection.append(artifact3)
```

While collections are indirectly annotated through the annotations of the artifacts they contain, you can also add collection-level annotations. Like artifacts, collections link to projects, runs, ulabels, records, and most other registries.

[^starsnowflake]: You can consider the SQL table underlying {class}`~lamindb.Artifact` your _fact table_ and all other tables for other entities your _dimension tables_ in a star or Snowflake schema ([see Wikipedia](https://en.wikipedia.org/wiki/Fact_table)).

[^protectproject]: The project annotation of the artifact is protected against the deletion of the project. If a user with necessary rights attempts to delete the project, they will get an error.
