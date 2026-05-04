# Organize datasets

```{toctree}
:maxdepth: 1
:hidden:

curate
```

This guide walks through organizing datasets with files & folders, with annotations in a database, and with versioned collections.

## Organize datasets as files and folders

If a database seems daunting, you can think of and use lamindb like a versioned file system in which you organize artifacts into virtual folders by using `/`-seperated keys, similar to AWS S3. For a single file, you'd call:

```python
artifact1 = ln.Artifact("./dataset.csv", key="project1/dataset1.csv").save()  # ingest the file in "folder" project/
```

For convenience, if you want to create an artifact for every file in a directory, use :meth:`~lamindb.Artifact.from_dir`:

```python
artifacts = ln.Artifact.from_dir("./project1/").save()  # create one artifact per file in the directory
```

You can then query for all artifacts in the folder `"./project1/"` folder via:

```python
artifacts = ln.Artifact.filter(key__startswith="project1/")  # query artifacts via the folder prefix
```

:::{dropdown} What if I do not care about the metadata and version of every file in a folder?

In some cases a folder _is_ the dataset and you do not care about fine-grained information for every file in it.
Consider then saving the whole directory as a single artifact:

```python
ln.Artifact("./folder_abc", key="folder_abc").save()  # create a single artifact for the whole "folder_abc/" directory
```

:::

## Organize datasets with annotations in a database

You can consider the {class}`~lamindb.Artifact` registry your central registry and all other registries as context that allow you to find & query artifacts based on different entities you care about.[^starsnowflake]

### Auto-generated annotations

Several {class}`~lamindb.Artifact` fields like `created_by`, `created_at`, `size`, etc. are automatically populated and you can use them to retrieve a set of artifacts:

```python
artifacts = ln.Artifact.filter(
    created_by__handle="falexwolf",  # created by user with handle falexwolf
    created_at__gt="2023-06-24",  # created after June 24th, 2023
    size__lt=1e9,   # smaller than 1GB
    description__icontains="scrnaseq",   # containing "scrnaseq" in description, case-insensitive
    suffix=".parquet",   # with a .parquet suffix
    n_observations__gt=1000,   # with more than 1000 observations
    n_files__gt=1000,    # folder-like artifacts with more than 1000 files
    otype="DataFrame",   # that are DataFrames
    created_on__name="my-branch",   # created on a specific branch or environment
    run=run,   # created by a specific run
    transform__name="my-script.py",   # created by a specific script/notebook
)
```

### Organizing artifacts with projects

What if an artifact is relevant to **multiple projects**?
A dataset that's in the `project1/` folder as in the example above cannot **also** be in a `project2/` folder.
You can solve this problem by annotating the artifact with projects:

```python
project1 = ln.Project(name="Project 1").save()  # create project 1
project2 = ln.Project(name="Project 2").save()  # create project 2
artifact1.projects.add(project1, project2)  # annotate dataset1
```

This allows you to find `artifact1` when you're querying for the datasets in any project that labels it.
For example, `artifact1` will be in the query results of both of these queries:

```python
artifacts_in_project1 = ln.Artifact.filter(projects=project1)  # all datasets in project1
artifacts_in_project2 = ln.Artifact.filter(projects=project2)  # all datasets in project2
```

Another advantage of this is that you don't have trust file paths anymore.
A folder structure in a file path might be renamed, and then your retrieval logic breaks.
A project that annotates an artifact or another object **cannot** be deleted[^fkprotect] so you can always trust that the query succeeds.

## Organizing artifacts with other labels

Often times you also want to annotate with other entities, not just projects. LaminDB offers two main classes for this: {class}`~lamindb.Record` for metadata records and {class}`~lamindb.ULabel` for simple labels. You can use these together with entities in modules such as {mod}`bionty` in full analogy with `Project`. For example:

```python
import bionty as bt

ulabel1 = ln.ULabel(name="raw_data").save()
record1 = ln.Record(name="My sample", features={"gc_content": 0.5}).save()
cell_type1 = bt.CellType.from_source(name="T cell").save()
artifact1.ulabels.add(ulabel1)
artifact1.records.add(record1)
artifact1.cell_types.add(cell_type1)
```

### Annotating with features

To annotate with non-categorical data types or to disambiguate categorical annotations, you might want consider annotations with features, using {class}`~lamindb.Feature` objects.

```python
experiment_type = ln.Record.get(name="Experiments")
ln.Feature(name="gc_content", dtype=float).save()
ln.Feature(name="experiment", dtype=experiment_type).save()
```

During annotation, feature names and data types are validated against these definitions.

```python
artifact.features.set_values({
    "gc_content": 0.55,
    "experiment": "Experiment 1",  # needs to exist under the "Experiments" record type
})
```

### Auto-annotating based on parsed metadata

When you work with structured data formats like `DataFrame`, `AnnData`, or similar, it often makes sense to validate their content. During validation, the content is parsed and can hence be used for annotation. This behavior is triggered if you pass a {class}`~lamindb.Schema` to {class}`~lamindb.Artifact`.

```python
# validate columns in the dataframe and map them on features, auto-annotate with parse metadata
ln.Artifact.from_dataframe(df, schema="valid_features").save()
```

## Publishing versioned collections of artifacts

In some cases, you don't just want to group a set of artifacts by different dimensions of metadata, but you also want to version the set. For this, you can use {class}`~lamindb.Collection`:

```python
collection = ln.Collection([artifact1, artifact2], key="my_dataset_release").save()
```

Unlike sets of artifacts defined through folders or through metadata annotations, where it's possible that artifacts are added or removed, a `collection` _guarantees_ that you get one exact immutable set of artifacts.

Artifacts are versiond based on the hash of their content. Collections are versioned based on the top-level hash of its artifact hashes.

[^fkprotect]: Its foreign key is protected in the `ArtifactProject` link model.

[^starsnowflake]: You can consider the SQL table underlying {class}`~lamindb.Artifact` your _fact table_ and all other tables for other entities your _dimension tables_ in a star or Snowflake schema ([see Wikipedia](https://en.wikipedia.org/wiki/Fact_table)).
