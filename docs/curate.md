# Curate datasets

Data curation with LaminDB ensures your datasets are **validated** and **queryable**. This guide shows you how to transform data into clean, annotated datasets.

Curating a dataset with LaminDB means three things:

- **Validate** that the dataset matches a desired schema.
- **Standardize** the dataset (e.g., by fixing typos, mapping synonyms) or update registries if validation fails.
- **Annotate** the dataset by linking it against metadata entities so that it becomes queryable.

In this guide we'll curate common data structures. Here is a [guide](/faq/curate-any) for the underlying low-level API.

Note: If you know either `pydantic` or `pandera`, here is an [FAQ](/faq/pydantic-pandera) that compares LaminDB with both of these tools.

```python
# pip install lamindb
!lamin init --storage ./test-curate --modules bionty
```

```python
import lamindb as ln

ln.track()
```

<!-- #region -->

## Schema design patterns

A {class}`~lamindb.Schema` in LaminDB is a specification that defines the expected structure, data types, and validation rules for a dataset.
It is similar to `pydantic.Model` for dictionaries, and `pandera.Schema`, and `pyarrow.lib.Schema` for tables, but supporting more complicated data structures.

Schemas ensure data consistency by defining:

- What {class}`~lamindb.Feature`s (dimensions) exist in your dataset
- What data types those features should have
- What values are valid for categorical features
- Which {class}`~lamindb.Feature`s are required vs optional

An exemplary schema:

```python
schema = ln.Schema(
    name="experiment_schema",           # human-readable name
    features=[                          # required features
        ln.Feature(name="cell_type", dtype=bt.CellType),
        ln.Feature(name="treatment", dtype=str),
    ],
    otype="DataFrame"                   # object type (DataFrame, AnnData, etc.)
)
```

For composite data structures using slots:

```{dropdown} What are slots?

For composite data structures, you need to specify which component contains which schema, for example, to validate both cell metadata in `.obs` and gene metadata in `.var` within the same schema.
Each slot is a key like `"obs"` for AnnData observations,`"rna:var"` for MuData modalities, or `"attrs:nested:key"` for SpatialData annotations.
```

```python
# AnnData with multiple "slots"
adata_schema = ln.Schema(
    otype="AnnData",
    slots={
        "obs": cell_metadata_schema,     # cell annotations
        "var.T": gene_id_schema          # gene-derived features
    }
)
```

Before diving into curation, let's understand the different schema approaches and when to use each one.
Think of schemas as rules that define what valid data should look like.

<!-- #endregion -->

### Flexible schema

Use when: You want to validate those columns whose names match feature names in your `Feature` registry.

```{eval-rst}
.. literalinclude:: scripts/define_valid_features.py
   :language: python
```

### Minimal required schema

Use when: You need certain columns but want flexibility for additional metadata.

```{eval-rst}
.. literalinclude:: scripts/define_mini_immuno_schema_flexible.py
   :language: python
```

<!-- #region -->

### Strict Schema

Use when: You need complete control over data structure and values.

```python
# Only allows specified columns
schema = ln.Schema(
    features=[...],
    minimal_set=True,  # whether all passed features are required
    maximal_set=False  # whether additional features are allowed
)
```

<!-- #endregion -->

## DataFrame

### Step 1: Load and examine your data

We'll be working with the mini immuno dataset:

```python
df = ln.examples.datasets.mini_immuno.get_dataset1(
    with_cell_type_synonym=True, with_cell_type_typo=True
)
df
```

### Step 2: Set up your metadata registries

Before creating a schema, ensure your registries have the right features and labels:

```{eval-rst}
.. literalinclude:: scripts/define_mini_immuno_features_labels.py
   :language: python
```

### Step 3: Create your schema

```python
schema = ln.examples.datasets.mini_immuno.define_mini_immuno_schema_flexible()
schema.describe()
```

<!-- #region -->

### Step 4: Initialize Curator and first validation

If you expect the validation to pass, you can directly register an artifact by providing the schema:

```python

artifact = ln.Artifact.from_dataframe(df, key="examples/my_curated_dataset.parquet", schema=schema).save()
```

<!-- #endregion -->

The {meth}`~lamindb.curators.core.Curator.validate` method validates that your dataset adheres to the criteria defined by the `schema`.
It identifies which values are already validated (exist in the registries) and which are potentially problematic (do not yet exist in our registries).

```python
try:
    curator = ln.curators.DataFrameCurator(df, schema)
    curator.validate()
except ln.errors.ValidationError as error:
    print(error)
```

### Step 5: Fix validation issues

```python
# check the non-validated terms
curator.cat.non_validated
```

For `cell_type_by_expert`, we saw 2 terms are not validated.

First, let's standardize synonym "B-cell" as suggested

```python
curator.cat.standardize("cell_type_by_expert")
```

```python
# now we have only one non-validated cell type left
curator.cat.non_validated
```

For "CD8-pos alpha-beta T cell", let's understand which cell type in the public ontology might be the actual match.

```python
# to check the correct spelling of categories, pass `public=True` to get a lookup object from public ontologies
# use `lookup = curator.cat.lookup()` to get a lookup object of existing records in your instance
lookup = curator.cat.lookup(public=True)
lookup
```

```python
# here is an example for the "cell_type" column
cell_types = lookup["cell_type_by_expert"]
cell_types.cd8_positive_alpha_beta_t_cell
```

```python
# fix the cell type name
df["cell_type_by_expert"] = df["cell_type_by_expert"].cat.rename_categories(
    {"CD8-pos alpha-beta T cell": cell_types.cd8_positive_alpha_beta_t_cell.name}
)
```

For perturbation, we want to add the new values: "DMSO", "IFNG"

```python
# this adds perturbations that were _not_ validated
curator.cat.add_new_from("perturbation")
```

```python
ln.Feature.get(name="perturbation")
```

```python
# validate again
curator.validate()
```

### Step 6: Save your curated dataset

```python
artifact = curator.save_artifact(key="examples/my_curated_dataset.parquet")
```

```python
artifact.describe()
```

## Common fixes

This section covers the most frequent curation issues and their solutions.
Use this as a reference when validation fails.

### Feature validation issues

<!-- #region -->

**Issue**: "Column not in dataframe"

```
"column 'treatment' not in dataframe. Columns in dataframe: ['drug', 'timepoint', ...]"
```

**Solutions**:

```python
# Solution 1: Rename columns to match schema
df = df.rename(columns={
    'treatment': 'drug',
    'time': 'timepoint',
    ...
})

# Solution 2: Create missing columns
df['treatment'] = 'unknown'  # Add with default value (or define Feature.default_value)

# Solution 3: Modify schema to match your data
schema = ln.Schema(
    features=[
        ln.Feature.get(name="drug"),  # Use actual column name
        ln.Feature.get(name="timepoint"),
    ],
    ...
)
```

<!-- #endregion -->

### Value validation issues

<!-- #region -->

**Issue**: "Terms not validated in feature 'perturbation'"

```
2 terms not validated in feature 'cell_type': 'B-cell', 'CD8-pos alpha-beta T cell'
    1 synonym found: "B-cell" → "B cell"
    → curate synonyms via: .standardize("cell_type")
    for remaining terms:
    → fix typos, remove non-existent values, or save terms via: curator.cat.add_new_from('cell_type')
```

**Solutions**:

```python
# Solution 1: Use automatic standardization if given hint (handles synonyms))
curator.cat.standardize('cell_type')

# Solution 2: Manual mapping for complex cases
value_mapping = {
    'T-cells': 'T cell',
    'B-cells': 'B cell',
}
df['cell_type'] = df['cell_type'].map(value_mapping).fillna(df['cell_type'])

# Solution 3: Use public ontology lookup for correct names
lookup = curator.cat.lookup(public=True)
cell_types = lookup["cell_type"]
df['cell_type'] = df['cell_type'].cat.rename_categories({
    'CD8-pos T cell': cell_types.cd8_positive_alpha_beta_t_cell.name
})

# Solution 4: Add new legitimate terms
curator.cat.add_new_from("cell_type")
```

<!-- #endregion -->

### Data type issues

<!-- #region -->

<cell_type>markdown</cell_type>**Issue**: "Expected categorical data, got object"

```
TypeError: Expected categorical data for cell_type, got object
```

**Solutions**:

```python
# Solution 1: Convert to categorical
df['cell_type'] = df['cell_type'].astype('category')

# Solution 2: Use coercion in feature definition
ln.Feature(name="cell_type", dtype=bt.CellType, coerce=True).save()
```

<!-- #endregion -->

## External data validation

Since not all metadata is always stored within the dataset itself, it is also possible to validate external metadata.

```{eval-rst}
.. literalinclude:: scripts/curate_dataframe_external_features.py
   :language: python
   :caption: curate_dataframe_external_features.py
```

```python
!python scripts/curate_dataframe_external_features.py
```

## AnnData

`AnnData` like all other data structures that follow is a composite structure that stores different arrays in different `slots`.

### Allow a flexible schema

We can also allow a flexible schema for an `AnnData` and only require that it's indexed with Ensembl gene IDs.

```{eval-rst}
.. literalinclude:: scripts/curate_anndata_flexible.py
   :language: python
   :caption: curate_anndata_flexible.py
```

Let's run the script.

```python
!python scripts/curate_anndata_flexible.py
```

Under-the-hood, this uses the following build-in schema ({func}`~lamindb.examples.schemas.anndata_ensembl_gene_ids_and_valid_features_in_obs`):

```{eval-rst}
.. literalinclude:: scripts/define_schema_anndata_ensembl_gene_ids_and_valid_features_in_obs.py
   :language: python
```

This schema tranposes the `var` DataFrame during curation, so that one validates and annotates the columns of `var.T`, i.e., `[ENSG00000153563, ENSG00000010610, ENSG00000170458]`.
If one doesn't transpose, one would annotate the columns of `var`, i.e., `[gene_symbol, gene_type]`.

```{eval-rst}
.. image:: https://lamin-site-assets.s3.amazonaws.com/.lamindb/gLyfToATM7WUzkWW0001.png
    :width: 800px
```

### Fix validation issues

```python
adata = ln.examples.datasets.mini_immuno.get_dataset1(
    with_gene_typo=True, with_cell_type_typo=True, otype="AnnData"
)
adata
```

```python
schema = ln.examples.schemas.anndata_ensembl_gene_ids_and_valid_features_in_obs()
schema.describe()
```

Check the slots of a schema:

```python
schema.slots
```

```python
curator = ln.curators.AnnDataCurator(adata, schema)
try:
    curator.validate()
except ln.errors.ValidationError as error:
    print(error)
```

As above, we leverage a lookup object with valid cell types to find the correct name.

```python
valid_cell_types = curator.slots["obs"].cat.lookup()["cell_type_by_expert"]
adata.obs["cell_type_by_expert"] = adata.obs[
    "cell_type_by_expert"
].cat.rename_categories(
    {"CD8-pos alpha-beta T cell": valid_cell_types.cd8_positive_alpha_beta_t_cell.name}
)
```

The validated `AnnData` can be subsequently saved as an {class}`~lamindb.Artifact`:

```python
adata.obs.columns
```

```python
curator.slots["var.T"].cat.add_new_from("columns")
```

```python
curator.validate()
```

```python
artifact = curator.save_artifact(key="examples/my_curated_anndata.h5ad")
```

Access the schema for each slot:

```python
artifact.features.slots
```

The saved artifact has been annotated with validated features and labels:

```python
artifact.describe()
```

## Unstructured dictionaries

Most datastructures support unstructured metadata stored as dictionaries:

- Pandas DataFrames: `.attrs`
- AnnData: `.uns`
- MuData: `.uns` and `modality:uns`
- SpatialData: `.attrs`

Here, we exemplary show how to curate such metadata for AnnData:

```{eval-rst}
.. literalinclude:: scripts/define_schema_anndata_uns.py
   :language: python
   :caption: define_schema_anndata_uns.py
```

```python
!python scripts/define_schema_anndata_uns.py
```

```{eval-rst}
.. literalinclude:: scripts/curate_anndata_uns.py
   :language: python
   :caption: curate_anndata_uns.py
```

```python
!python scripts/curate_anndata_uns.py
```

## MuData

```{eval-rst}
.. literalinclude:: scripts/curate_mudata.py
   :language: python
   :caption: curate_mudata.py
```

```python
!python scripts/curate_mudata.py
```

## SpatialData

```{eval-rst}
.. literalinclude:: scripts/define_schema_spatialdata.py
   :language: python
   :caption: define_schema_spatialdata.py
```

```python
!python scripts/define_schema_spatialdata.py
```

```{eval-rst}
.. literalinclude:: scripts/curate_spatialdata.py
   :language: python
   :caption: curate_spatialdata.py
```

```python
!python scripts/curate_spatialdata.py
```

## TiledbsomaExperiment

```{eval-rst}
.. literalinclude:: scripts/curate_soma_experiment.py
   :language: python
   :caption: curate_soma_experiment.py
```

```python
!python scripts/curate_soma_experiment.py
```

## Other data structures

If you have other data structures, read: {doc}`/faq/curate-any`.

```python
!rm -rf ./test-curate
!rm -rf ./small_dataset.tiledbsoma
!lamin delete --force test-curate
```
