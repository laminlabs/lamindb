# Basic organization

It is possible to simply treat `LaminDB` as a file system and be content with organizing datasets as artifacts in which `key` represents virtual folders:

```python
ln.Artifact("my_dataset.csv", key="my_project/my_dataset.csv).save()
```

Here is a recommended workflow for ensuring datasets in LaminDB are FAIR.

| Step                 | Action                                               | Primitive                                                                                                           |
| :------------------- | :--------------------------------------------------- | :------------------------------------------------------------------------------------------------------------------ |
| **1. Define Scope**  | Group work by study or budget                        | {class}`Project`                                                                                                    |
| **2. Register Data** | Track individual files/arrays and lineage            | {class}`Artifact`                                                                                                   |
| **3. Annotate**      | Annotate with features, metadata records, and labels | {class}`Feature`, {class}`Record` {class}`ULabel`, {class}`bionty.Disease`, {class}`bionty.ExperimentalFactor`, ... |
| **4. Publish**       | Release an immutable, versioned bundle               | {class}`Collection`                                                                                                 |

Transitioning from folder-based structures to LaminDB involves a shift from **physical organization** (file paths) to **logical organization** (metadata and registries). While storage locations can still be managed via key prefixes, the registry acts as the primary interface for data access and discovery.

### Physical vs. Logical Organization

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
