# Group, annotate, and version datasets

This guide walks from using LaminDB like a file system to using it as a publishing system for curated data products.

## Using `lamindb` like a file system

It is possible to treat `lamindb` as a file system and and organize datasets as artifacts in which a `/`-seperated `key` represents virtual folders, just like on AWS S3. For a single file, you'd call:

```python
ln.Artifact("my_dataset.csv", key="my_project/my_dataset.csv").save()
```

To create an artifact for every file in a directory, use :meth:`~lamindb.Artifact.from_dir`:

```python
artifacts = ln.Artifact.from_dir("project_alpha/run_001").save()  # create one artifact per file in the directory
artifacts = ln.Artifact.filter(key__startswith="project_alpha/run_001/")  # query ingested artifacts via the folder prefix
```

In some cases, you might not care about metadata and lineage of a file in a directory. In this case, you can also save a whole directory as a single artifact:

```python
ln.Artifact("project_alpha/run_001", key="project_alpha/run_001").save()  # create a single artifact for the whole directory
```

## Basic annotations of datasets

Beyond a `key` or storage path, the most basic way to organize an artifact is by providing a `description` and a `version`:

```python
artifact = ln.Artifact(
    "my_dataset.csv",
    description="A dataset from run 001",
    version="1.0",
).save()
```

You can further annotate datasets with labels. For instance, using the built-in label registry {class}`~lamindb.ULabel`:

```python
label = ln.ULabel(name="Project Alpha").save()
artifact.ulabels.add(label)
```

This allows you to easily discover datasets later without relying on their file paths:

```python
artifacts = ln.Artifact.filter(ulabels__name="Project Alpha").all()
```

For comprehensive schema validation and standardizing metadata using biological ontologies (e.g., cell types, genes), see the [curate guide](curate.md).

## Publishing versioned data products

Here is a recommended workflow for ensuring datasets in LaminDB are FAIR.

| Step                 | Action                                               | Primitive                                                                                                           |
| :------------------- | :--------------------------------------------------- | :------------------------------------------------------------------------------------------------------------------ |
| **1. Define Scope**  | Group work by study or budget                        | {class}`Project`                                                                                                    |
| **2. Register Data** | Track individual files/arrays and lineage            | {class}`Artifact`                                                                                                   |
| **3. Annotate**      | Annotate with features, metadata records, and labels | {class}`Feature`, {class}`Record` {class}`ULabel`, {class}`bionty.Disease`, {class}`bionty.ExperimentalFactor`, ... |
| **4. Publish**       | Release an immutable, versioned bundle               | {class}`Collection`                                                                                                 |

Transitioning from folder-based structures to LaminDB involves a shift from **physical organization** (file paths) to **logical organization** (metadata context). While storage locations can still be managed via key prefixes, the `Artifact` and its related registries act as the primary interface for data access and discovery.

### Physical vs. logical organization

In traditional file systems, context is embedded in nested folder hierarchies. In LaminDB, data is organized via three core primitives: `Artifact`, `Collection`, and `Project`.

#### 1. Artifacts (Atomic Data)

An `Artifact` represents a discrete data object, such as a `.zarr` store, `.csv` file, or `.pdf` report.

- **Mapping:** Replaces individual files or sub-folders.
- **Property:** Each artifact is versioned and associated with biological or technical metadata (e.g., cell type, gene, or assay). The unique ID and metadata are the primary identifiers; the physical S3 path is secondary.

#### 2. Collections (Data Products)

A `Collection` is a versioned, immutable group of `Artifacts`.

- **Mapping:** Replaces "release" or "publication" folders.
- **Use Case:** Use collections to publish curated datasets. Unlike folders, collections are stable; a specific version of a collection always references the same set of artifact versions, ensuring reproducibility.

#### 3. Projects (Administrative Scope)

A `Project` provides an administrative boundary for research efforts.

- **Mapping:** Replaces top-level study or department directories.
- **Function:** Tracks the lineage of `Transforms` (notebooks/pipelines), `Artifacts`, and `Collections` created for a specific research goal.

### Storage Paths and Key Prefixes

LaminDB can track data where it lives (external) or manage it within a storage location. While you can specify **key prefixes** (virtual folders) to maintain compatibility with legacy scripts that require specific S3 paths, organization should happen within the registry. Accessing data via `ln.Artifact.filter()` or `ln.Collection.get()` is preferred over hardcoding file paths.
