---
execute_via: python
---

# Tutorial

This tutorial complements the [quickstart](https://docs.lamin.ai/introduction#quickstart): there, you skim through core features; here, you learn the core operations end-to-end with conceptual context.

## Track changes

If you don't have a LaminDB instance, create one in your current working directory on the terminal:

<!-- #skip_laminr -->

```bash
lamin init --modules bionty
```

<!-- #end_skip_laminr -->

<!-- #region -->

:::{dropdown} Via the R shell

```R
library(laminr)
lc <- import_module("lamin_cli")
lc$init(storage = "./lamin-tutorial", modules = "bionty")
```

:::

:::{dropdown} What else can I configure during setup?

1. You can pass a cloud storage location to `--storage` (S3, GCP, R2, HF, etc.)
   ```python
   --storage s3://my-bucket
   ```
2. Instead of the default SQLite database pass a Postgres connection string to `--db`:
   ```python
   --db postgresql://<user>:<pwd>@<hostname>:<port>/<dbname>
   ```
3. Instead of a default instance name derived from the storage location or the Postgres database, provide a custom name:
   ```python
   --name my-name
   ```
4. Mount additional modules:
   ```python
   --modules bionty,pertdb,custom1
   ```

For more info: {doc}`/setup`

:::

<!-- #endregion -->

### Track a notebook or script

Let's now track the notebook or script that is running.

```python
import lamindb as ln

ln.track()  # track the current notebook or script
```

Calling {meth}`~lamindb.track` links this execution to the data you'll save in this run. You can inspect your scripts & notebooks in the {class}`~lamindb.Transform` registry and their runs in the {class}`~lamindb.Run` registry.

```python
ln.Transform.to_dataframe()
```

```python
ln.Run.to_dataframe()
```

:::{dropdown} What happened under the hood?

1. The full run environment and imported package versions of current notebook were detected
2. Notebook metadata was detected and stored in a {class}`~lamindb.Transform` record with a unique id
3. Run metadata was detected and stored in a {class}`~lamindb.Run` record with a unique id

The {class}`~lamindb.Transform` registry stores data transformations: scripts, notebooks, pipelines, functions.

The {class}`~lamindb.Run` registry stores executions of transforms. Many runs can be linked to the same transform if executed with different context (time, user, input data, etc.).

For more info: {doc}`/track`

:::

<!-- #region -->

:::{dropdown} How do I track a workflow or a pipeline instead of a notebook or script?

Use {meth}`~lamdindb.flow`:

```python
import lamindb as ln


@ln.flow()
def ingest_dataset(key: str) -> ln.Artifact:
    df = ln.examples.datasets.mini_immuno.get_dataset1()
    return ln.Artifact.from_dataframe(df, key=key).save()


ingest_dataset(key="my_analysis/dataset.parquet")
```

For more info: {doc}`/track`

Leverage a pipeline integration, see: {doc}`docs:pipelines`.

You can also manually create a transform and manage its runs with your custom logic:

```python
transform = ln.Transform(key="my_pipeline", version="1.2.0")
```

:::

<!-- #endregion -->

:::{dropdown} Why should I care about emergent data lineage?

When humans & agents are in the loop, most mistakes happen outside deterministic computational pipelines. They happen in between pipelines and when running notebooks & scripts.

{func}`~lamindb.track` makes notebooks & derived results reproducible & auditable, enabling to learn from mistakes.

:::

:::{dropdown} Is this compliant with OpenLineage?

Yes. What OpenLineage calls a "job", LaminDB calls a "transform". What OpenLineage calls a "run", LaminDB calls a "run".

:::

### Manage changes with branches

When you need more control you can also use git-inspired change management via contribution branches that you merge into main.

<!-- #region -->

```bash
lamin switch -c my_branch
# ... make changes ...
lamin switch main
lamin merge my_branch
```

<!-- #endregion -->

For more info: {class}`~lamindb.Branch`

### Manage agent plans

If you're working with an agent plan, here is a way to save it:

<!-- #region -->

```bash
lamin save /path/to/.cursor/plans/my-agent-plan.md
```

It you switched to a branch, the plan will be linked to that branch. Here is how to link a plan to a run:

```python
ln.track(plan=".plans/my-agent-plan.md")  # link plan artifact to this run
```

<!-- #endregion -->

## Manage datasets

{class}`~lamindb.Artifact` is the core object for datasets and models in LaminDB, whether they are files, folders, tables, or array-like data.

It gives you one interface for registration, access, annotation, lineage, and queries across local storage, AWS S3 (`s3://`), Google Cloud (`gs://`), Hugging Face (`hf://`), and any other file system supported by `fsspec`.

### Create an artifact

Let's first look at an exemplary `DataFrame`. We'll cover simple files, folders, and the `AnnData` format later in this tutorial.

```python
df = ln.examples.datasets.mini_immuno.get_dataset1(with_typo=True)
df
```

:::

This is how you create an artifact from a `DataFrame`.

```python
artifact = ln.Artifact.from_dataframe(df, key="my_datasets/rnaseq1.parquet").save()
artifact.describe()
```

:::{dropdown} Which fields are populated when creating an artifact?

Basic fields:

- `uid`: universal ID
- `key`: a (virtual) relative path of the artifact in `storage`
- `description`: an optional string description
- `storage`: the storage location (the root, say, an S3 bucket or a local directory)
- `suffix`: an optional file/path suffix
- `size`: the artifact size in bytes
- `hash`: a hash useful to check for integrity and collisions (is this artifact already stored?)
- `created_at`: time of creation
- `updated_at`: time of last update

Provenance-related fields:

- `created_by`: the {class}`~lamindb.User` who created the artifact
- `run`: the {class}`~lamindb.Run` of the {class}`~lamindb.Transform` that created the artifact

For a full reference, see {class}`~lamindb.Artifact`.

:::

:::{dropdown} How does LaminDB compare to a AWS S3?

LaminDB provides a database on top of AWS S3 (or GCP storage, file systems, etc.).

Similar to organizing files with paths, you can organize artifacts using the `key` parameter of {class}`~lamindb.Artifact`.

However, you'll see that you can more conveniently query data by entities you care about: people, code, experiments, genes, proteins, cell types, etc.

:::

:::{dropdown} What exactly happens during save?

In the database: A SQL record is inserted into the `Artifact` registry. If the SQL record exists already based on comparing its `hash`, it's returned.

In storage:

- If the default storage is in the cloud, `.save()` triggers an upload for a local artifact.
- If the artifact is in a registered storage location, only the metadata of the artifact is saved as a SQL record to the `Artifact` registry.

:::

### Access artifacts

Get the artifact by `key`:

```python
artifact = ln.Artifact.get(key="my_datasets/rnaseq1.parquet")
```

:::

Load it into memory:

```python
artifact.load()
```

Typically your artifact is in a cloud storage location. To get a local file path to it, call {meth}`~lamindb.Artifact.cache`:

```python
artifact.cache()
```

If the data is large, you might not want to cache but stream it via {meth}`~lamindb.Artifact.open`. For more: {doc}`arrays`

### Update & delete artifacts

Here is how to update & delete an artifact:

```python
artifact.description = "My updated description"
artifact.save()  # persist metadata changes

artifact.delete()  # move to trash
artifact.restore()  # restore from trash

# artifact.delete(permanent=True)  # permanently delete
```

:::{dropdown} What happens when I delete an artifact?

By default, deleting moves an artifact into the trash so it no longer appears in standard queries.

To restore it, call `artifact.restore()`. For archive/trash semantics and query patterns, see {doc}`manage-changes`.

:::

### Trace data lineage

You can understand where an artifact comes from by looking at its {class}`~lamindb.Transform` & {class}`~lamindb.Run` attributes:

```python
artifact.transform
```

```python
artifact.run
```

Or visualize deeper data lineage with the `view_lineage()` method. Here we're only one step deep:

```python
artifact.view_lineage()
```

:::::{dropdown} Show me a more interesting example, please!

::::{tab-set}
:::{tab-item} Py

Explore and load the notebook [here](https://lamin.ai/laminlabs/lamindata/transform/F4L3oC6QsZvQ).

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/KQmzmmLOeBN0C8Yk0004.png" width="800">
:::
:::{tab-item} Hub

Explore data lineage interactively [here](https://lamin.ai/laminlabs/lamindata/artifact/W1AiST5wLrbNEyVq).

```{raw} html
<video width="500" controls>
  <source src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/Xdiikc2c1tPtHcvF0000.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>
```

:::
::::
:::::

Once you're done, at the end of your notebook or script, call {meth}`~lamindb.finish`. Here, we're not yet done so we're commenting it out.

```python
# ln.finish()  # mark run as finished, save execution report, source code & environment
```

:::{dropdown} Saving reports when using R

If you did _not_ use RStudio's notebook mode, you have to render an HTML externally.

1. Render the notebook to HTML via one of:

   - In RStudio, click the "Knit" button
   - From the command line, run

     ```bash
     Rscript -e 'rmarkdown::render("tutorial.Rmd")'
     ```

   - Use the `rmarkdown` package in R

     ```r
     rmarkdown::render("tutorial.Rmd")
     ```

2. Save it via one of:

   - the `save()` command in the `lamin_cli` module from R

     ```r
     lc <- import_module("lamin_cli")
     lc$save("tutorial.Rmd")
     ```

   - the `lamin` CLI

     ```bash
     lamin save tutorial.Rmd
     ```

:::

<!-- #region -->

To create a new version of a notebook or script, run `lamin load` on the terminal, e.g.,

```bash
$ lamin load https://lamin.ai/laminlabs/lamindata/transform/13VINnFk89PE0004
→ notebook is here: mcfarland_2020_preparation.ipynb
```

<!-- #endregion -->

Note that data lineage also helps to understand what a dataset is being used for, not only where it comes from. Many datasets are being used over and over for different purposes and it's often useful to understand how.

### Annotate an artifact

You can annotate artifacts with features and labels. Features are measurement dimensions (e.g. `"organism"`, `"temperature"`) and labels are measured categories (e.g. `"human"`, `"mouse"`).

Let's annotate an artifact with a {class}`~lamindb.ULabel`, a built-in universal label registry.

```python
# create & save a universal label
my_experiment = ln.ULabel(name="My experiment").save()

# annotate the artifact with a universal label
artifact.ulabels.add(my_experiment)

# describe the artifact
artifact.describe()
```

:::

This is how you query artifacts based on the annotation.

```python
ln.Artifact.filter(ulabels=my_experiment).to_dataframe()
```

You can also annotate with labels from other registries, e.g., the biological ontologies in {mod}`bionty`.

```python
import bionty as bt

# create a cell type label from the source ontology
cell_type = bt.CellType.from_source(name="effector T cell").save()

# annotate the artifact with a cell type
artifact.cell_types.add(cell_type)

# describe the artifact
artifact.describe()
```

This is how you query artifacts by cell type annotations.

```python
ln.Artifact.filter(cell_types=cell_type).to_dataframe()
```

Here is how to update or remove annotations:

```python
# add labels
artifact.ulabels.add(my_experiment)
artifact.cell_types.add(cell_type)

# remove labels
artifact.ulabels.remove(my_experiment)
artifact.cell_types.remove(cell_type)
```

If you want to annotate by non-categorical metadata or indicate the feature for a label, annotate via features.

```python
# define the "temperature" & "experiment" features
temperature = ln.Feature(name="temperature", dtype=float).save()
experiment = ln.Feature(name="experiment", dtype=ln.ULabel).save()

# annotate the artifact
artifact.features.set_values({temperature: 21.6, experiment: "My experiment"})

# describe the artifact
artifact.describe()
```

This is how you query artifacts by features.

```python
ln.Artifact.filter(temperature=21.6).to_dataframe()
```

Optionally, after querying, this is how you clean up feature values from the artifact and remove the feature definitions.

```python
# remove feature values from this artifact
artifact.features.remove_values(["temperature", "experiment"])

# optionally remove feature definitions from the registry
ln.Feature.get(name="temperature").delete()
ln.Feature.get(name="experiment").delete()
```

### Validate an artifact

Validated datasets are more re-usable by analysts and machine learning models. You can define what a valid artifact should look like by defining a schema. If you do that, LaminDB will also automatically annotate with validated metadata which helps with finding the artifact.

:::{dropdown} Give me examples for findability and usability.

1. Findability: Which datasets measured expression of cell marker `CD14`? Which characterized cell line `K562`? Which have a test & train split? Etc.
2. Usability: Are there typos in feature names? Are there typos in labels? Are types and units of features consistent? Etc.

:::

```python
import bionty as bt  # <-- use bionty to access registries with imported public ontologies

# define a few more valid labels
ln.ULabel(name="DMSO").save()
ln.ULabel(name="IFNG").save()

# define a few more valid features
ln.Feature(name="perturbation", dtype=ln.ULabel).save()
ln.Feature(name="cell_type_by_model", dtype=bt.CellType).save()
ln.Feature(name="cell_type_by_expert", dtype=bt.CellType).save()
ln.Feature(name="assay_oid", dtype=bt.ExperimentalFactor.ontology_id).save()
ln.Feature(name="donor", dtype=str, nullable=True).save()
ln.Feature(name="sample_note", dtype=str, nullable=True).save()
ln.Feature(name="concentration", dtype=str).save()
ln.Feature(name="treatment_time_h", dtype="num", coerce=True).save()

# define a schema that merely enforces a feature identifier type
schema = ln.Schema(itype=ln.Feature).save()
```

If you pass a `schema` object to the `Artifact` constructor, the artifact will be validated & annotated. Let's try this.

```python
try:
    artifact = ln.Artifact.from_dataframe(
        df, key="my_datasets/rnaseq1.parquet", schema=schema
    )
except ln.errors.ValidationError as error:
    print(str(error))
```

Because there is a typo in the `perturbation` column, validation fails. Let's fix it by making a new version.

### Make a new version of an artifact

```python
# fix the "IFNJ" typo
df["perturbation"] = df["perturbation"].cat.rename_categories({"IFNJ": "IFNG"})

# create a new version
artifact = ln.Artifact.from_dataframe(
    df, key="my_datasets/rnaseq1.parquet", schema=schema
).save()

# see the annotations
artifact.describe()

# see all versions of the artifact
artifact.versions.to_dataframe()
```

The content of the dataset is now validated and the dataset is richly annotated and queryable by all entities that you defined.

<!-- #region -->

:::{dropdown} Can I also create new versions without passing `key`?

That works, too, you can use `revises`:

```python
artifact_v1 = ln.Artifact.from_dataframe(df, description="Just a description").save()
# below revises artifact_v1
artifact_v2 = ln.Artifact.from_dataframe(df_updated, revises=artifact_v1).save()
```

<br>

The good thing about passing `revises` is that you don't need to worry about coming up with naming conventions for paths. The good thing about versioning based on `key` is that it's how all data versioning tools are doing it and you can build a file hierarchy.

:::

<!-- #endregion -->

## Query & search registries

We've already seen a few queries. Let's now walk through the topic systematically.

To get an overview over all artifacts in your instance, call {class}`~lamindb.models.SQLRecord.df`.

```python
ln.Artifact.to_dataframe()
```

The `Artifact` registry additionally supports seeing all feature annotations of an artifact.

```python
ln.Artifact.to_dataframe(include="features")
```

LaminDB's central classes are registries that store {class}`~lamindb.models.SQLRecord` objects. If you want to see the fields of a registry, call `.describe()` on the class or auto-complete:

```python
ln.Artifact.describe()
```

Each registry is a table in the relational schema of the underlying database. With {func}`~lamindb.view`, you can see the latest additions to the database:

```python
ln.view()
```

:::{dropdown} Which registries have I already learned about?

- {class}`~lamindb.Artifact`: datasets & models stored as files, folders, or arrays
- {class}`~lamindb.Transform`: transforms of artifacts
- {class}`~lamindb.Run`: runs of transforms
- {class}`~lamindb.User`: users
- {class}`~lamindb.Storage`: local or cloud storage locations

:::

Every registry supports arbitrary relational queries using the class methods {class}`~lamindb.models.SQLRecord.get` and {class}`~lamindb.models.SQLRecord.filter`.
The syntax for it is Django's query syntax.

Here are some simple query examples.

```python
# get a single record (here the current notebook)
transform = ln.Transform.get(key="tutorial.ipynb")

# get a set of records by filtering for a directory (LaminDB treats directories like AWS S3, as the prefix of the storage key)
ln.Artifact.filter(key__startswith="my_datasets/").to_dataframe()

# query all artifacts ingested from a transform
artifacts = ln.Artifact.filter(transform=transform).all()

# query all artifacts ingested from a notebook with "tutor" in the description
artifacts = ln.Artifact.filter(
    transform__description__icontains="tutor",
).all()
```

:::{dropdown} What does a double underscore mean?

For any field, the double underscore defines a comparator, e.g.,

- `name__icontains="Martha"`: `name` contains `"Martha"` when ignoring case
- `name__startswith="Martha"`: `name` starts with `"Martha`
- `name__in=["Martha", "John"]`: `name` is `"John"` or `"Martha"`

For more info: {doc}`query-search`

:::

:::{dropdown} Can I chain filters and searches?

Yes: `ln.Artifact.filter(suffix=".jpg").search("my image")`

:::

The class methods {class}`~lamindb.models.SQLRecord.search` and {class}`~lamindb.models.SQLRecord.lookup` help with approximate matches.

```python
# search artifacts
ln.Artifact.search("iris").to_dataframe().head()

# search transforms
ln.Transform.search("tutor").to_dataframe()

# look up records with auto-complete
records = ln.Record.lookup()
```

:::{dropdown} Show me a screenshot

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/lgRNHNtMxjU0y8nIagt7.png" width="400px">

:::

For more info: {doc}`query-search`

## Manage folders & storage locations

Let's look at a folder in the cloud that contains 3 sub-folders storing images & metadata of Iris flowers, generated in 3 subsequent studies.

```python
# we use anon=True here in case no aws credentials are configured
ln.UPath("s3://lamindata/iris_studies", anon=True).view_tree()
```

Let's create an artifact for the first sub-folder.

```python
artifact = ln.Artifact("s3://lamindata/iris_studies/study0_raw_images").save()
artifact
```

As you see from {attr}`~lamindb.Artifact.path`, the folder was merely registered in its present storage location without copying it.

```python
artifact.path
```

LaminDB keeps track of all your storage locations.

```python
ln.Storage.to_dataframe()
```

<!-- #region -->

:::{dropdown} How do I create an artifact for a local file or folder?

Source path is local:

```python
ln.Artifact("./my_data.fcs", key="my_dataset.fcs")
ln.Artifact("./my_images/", key="my_images")
```

<br>

Upon `artifact.save()`, the source path will be copied or uploaded into your instance's current storage, visible & changeable via `ln.settings.storage`.

If the source path is remote _or_ already in a registered storage location (one that's registered in `ln.Storage`), `artifact.save()` will _not_ trigger a copy or upload but register the existing path.

```python
ln.Artifact("s3://my-bucket/my_dataset.fcs")  # key is auto-populated from S3, you can optionally pass a description
ln.Artifact("s3://my-bucket/my_images/")  # key is auto-populated from S3, you can optionally pass a description
```

<br>
You can use any storage location supported by `fsspec`.

:::

<!-- #endregion -->

<!-- #region -->

:::{dropdown} Are artifacts aware of array-like data?

Yes.

You can make artifacts from paths referencing array-like objects:

```python
ln.Artifact("./my_anndata.h5ad", key="my_anndata.h5ad")
ln.Artifact("./my_zarr_array/", key="my_zarr_array")
```

Or from in-memory objects:

```python
ln.Artifact.from_dataframe(df, key="my_dataframe.parquet")
ln.Artifact.from_anndata(adata, key="my_anndata.h5ad")
```

You can open large artifacts for slicing from the cloud or load small artifacts directly into memory via:

```python
artifact.open()
```

:::

<!-- #endregion -->

## Manage biological ontologies

Every {py:mod}`bionty` registry is based on configurable public ontologies (>20 of them) that are automatically leveraged during validation & annotation. Sometimes you want to access the public ontology directly.

```python
import bionty as bt

cell_type_ontology = bt.CellType.public()
cell_type_ontology
```

The returned object can be searched like you can search a registry.

```python
cell_type_ontology.search("gamma-delta T cell").head(2)
```

Because you can't update an external public ontology, you update the content of the corresponding registry. Here, you create a new cell type.

```python
# create an ontology-coupled cell type record and save it
neuron = bt.CellType.from_source(name="neuron").save()

# create a record to track a new cell state
new_cell_state = bt.CellType(
    name="my neuron cell state", description="explains X"
).save()

# express that it's a neuron state
new_cell_state.parents.add(neuron)

# view ontological hierarchy
new_cell_state.view_parents(distance=2)
```

For more info: {doc}`manage-ontologies`

## Manage AnnData artifacts

LaminDB supports a growing number of data structures: `DataFrame`, `AnnData`, `MuData`, `SpatialData`, and `Tiledbsoma` with their corresponding storage formats.

Let's go through the example of the quickstart, but store the dataset in an AnnData this time.

```python
# define var schema
var_schema = ln.Schema(itype=bt.Gene.ensembl_gene_id, dtype=int).save()

# define composite schema
anndata_schema = ln.Schema(
    otype="AnnData", slots={"obs": schema, "var.T": var_schema}
).save()
```

Validate & annotate an `AnnData`.

```python
import anndata as ad

# store the dataset as an AnnData object to distinguish data from metadata
adata = ad.AnnData(df.iloc[:, :3], obs=df.iloc[:, 3:-1])

# save curated artifact
artifact = ln.Artifact.from_anndata(
    adata, key="my_datasets/my_rnaseq1.h5ad", schema=anndata_schema
).save()
artifact.describe()
```

Because `AnnData` separates the high-dimensional count matrix that's typically indexed with Ensembl gene ids from the metadata, we're now working with two types of feature sets (`bt.Gene` for the counts and `ln.Feature` for the metadata). These correspond to the `obs` and the `var` schema in the `anndata_schema`.

If you want to find a dataset by whether it measured `CD8A`, you can do so as as follows.

```python
# query for all schemas that contain CD8A
schemas = ln.Schema.filter(genes__symbol="CD8A").all()

# query for all artifacts linked to these feature sets
ln.Artifact.filter(schemas__in=schemas).to_dataframe()
```

## Manage collections of datasets

How do you integrate new datasets with your existing datasets? Leverage {class}`~lamindb.Collection`.

```python
# a new dataset
df2 = ln.examples.datasets.mini_immuno.get_dataset2(otype="DataFrame")
adata = ad.AnnData(df2.iloc[:, :3], obs=df2.iloc[:, 3:-1])
artifact2 = ln.Artifact.from_anndata(
    adata, key="my_datasets/my_rnaseq2.h5ad", schema=anndata_schema
).save()
```

Create a collection using {class}`~lamindb.Collection`.

```python
collection = ln.Collection([artifact, artifact2], key="my-RNA-seq-collection").save()
collection.describe()
collection.view_lineage()
```

```python
# if it's small enough, you can load the entire collection into memory as if it was one
collection.load()

# typically, it's too big, hence, open it for streaming (if the backend allows it)
# collection.open()

# or iterate over its artifacts
collection.ordered_artifacts.all()

# or look at a DataFrame listing the artifacts
collection.ordered_artifacts.to_dataframe()
```

Directly train models on collections of `AnnData`.

<!-- #region -->
<!-- #skip_laminr -->

```python
# to train models, batch iterate through the collection as if it was one array
from torch.utils.data import DataLoader, WeightedRandomSampler
dataset = collection.mapped(obs_keys=["cell_medium"])
sampler = WeightedRandomSampler(
    weights=dataset.get_label_weights("cell_medium"), num_samples=len(dataset)
)
data_loader = DataLoader(dataset, batch_size=2, sampler=sampler)
for batch in data_loader:
    pass
```

<!-- #end_skip_laminr -->
<!-- #endregion -->

For more: {doc}`/arrays`
