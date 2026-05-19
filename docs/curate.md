---
execute_via: python
---

# Validate & standardize datasets

```{raw} html
<iframe width="560" height="315" src="https://www.youtube.com/embed/Ji6E7hTnReQ?si=K0OnU2MTGv4fIhFo" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
```

Curating a dataset means three things:

- **Validate** that the dataset matches a desired schema.
- If validation fails, **standardize** the dataset (e.g., by fixing typos, mapping synonyms) or update registries.
- **Annotate** the dataset by linking it against metadata entities so that it becomes queryable.

In other guides, we've mostly covered annotation. In this guide we'll curate common data structures focusing on validation and standardization.

<!-- #region -->

## Schema design patterns

A {class}`~lamindb.Schema` is a specification that defines the expected structure, data types, and validation rules for a dataset.
Unlike `pydantic.Model` for dictionaries or `pandera.Schema` / `pyarrow.lib.Schema` for tables, a `lamindb` schema enables queries in a database and supports complex data structures.
Here is an [FAQ](/faq/pydantic-pandera) that compares it with `pydantic` and `pandera`.

Schemas ensure data consistency by defining:

- Which features exist in your dataset
- What data types those features should have
- What values are valid for categorical features
- Which features are required vs optional

An exemplary schema that leverages {class}`~lamindb.Feature` objects to define features:

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

Or for a composite data structure, like an `AnnData`:

```python
schema = ln.Schema(
    otype="AnnData",
    slots={
        "obs": cell_metadata_schema,     # (pseudocode) cell annotations
        "var.T": gene_id_schema          # (pseudocode) gene-derived features
    }
)
```

```{dropdown} What are slots?

For composite data structures, you need to specify which component contains which schema, for example, to validate both cell metadata in `.obs` and gene metadata in `.var` within the same schema.
Each slot is a key like `"obs"` for AnnData observations,`"rna:var"` for MuData modalities, or `"attrs:nested:key"` for SpatialData annotations.
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
    features=[...],    # (pseudocode)
    minimal_set=True,  # whether all passed features are required
    maximal_set=False  # whether additional features are allowed
)
```

<!-- #endregion -->

## DataFrame

If you're not connected to a database, create one:

```python
!lamin init --storage ./test-curate --modules bionty
```

Let's import `lamindb` and optionally track this run:

```python
import lamindb as ln

ln.track()
```

We'll not be walking through the steps of curating a dataset.

### (1) Load and examine your dataset

We'll be working with the mini immuno dataset:

```python
df = ln.examples.datasets.mini_immuno.get_dataset1(
    with_cell_type_synonym=True, with_cell_type_typo=True
)
df
```

### (2) Set up your registries

Before creating a schema, ensure your registries have the right features and labels:

```{eval-rst}
.. literalinclude:: scripts/define_mini_immuno_features_labels.py
   :language: python
```

### (3) Create your schema

Let's instantiate the flexible schema we discussed earlier (available in our examples module):

```python
schema = ln.examples.datasets.mini_immuno.define_mini_immuno_schema_flexible()
schema.describe()
```

<!-- #region -->

### (4) Validate the dataset

:::{admonition} Shortcut
If you expect the validation to pass, you can directly ingest a validated artifact via:

```python
artifact = ln.Artifact.from_dataframe(df, key="examples/my_curated_dataset.parquet", schema=schema).save()
```

:::

<!-- #endregion -->

If you want full control over the validation process with access to standardization helpers, you can instantiate a `Curator` object. Its {meth}`~lamindb.curators.core.Curator.validate` method validates that your dataset adheres to the criteria defined by the `schema`.
It identifies which values are already validated (exist in the registries) and which are potentially problematic (do not yet exist in our registries).

```python
try:
    curator = ln.curators.DataFrameCurator(df, schema)
    curator.validate()
except ln.errors.ValidationError as error:
    print(error)
```

### (5) Fix validation errors

Check the non-validated terms:

```python
curator.cat.non_validated
```

For `cell_type_by_expert`, we see 2 terms are not validated.

First, let's standardize the synonym "B-cell" as suggested:

```python
curator.cat.standardize("cell_type_by_expert")
# now we have only one non-validated cell type left
curator.cat.non_validated
```

For "CD8-pos alpha-beta T cell", let's understand which cell type in the public ontology might be the actual match:

```python
# to check the correct spelling of categories, pass `public=True` to get a lookup object from public ontologies
# use `lookup = curator.cat.lookup()` to get a lookup object of existing records in your instance
lookup = curator.cat.lookup(public=True)
# here is an example for the "cell_type_by_expert" feature
cell_types = lookup["cell_type_by_expert"]
cell_types.cd8_positive_alpha_beta_t_cell
```

Now that we have a match, let's fix it:

```python
# fix the cell type name
df["cell_type_by_expert"] = df["cell_type_by_expert"].cat.rename_categories(
    {"CD8-pos alpha-beta T cell": cell_types.cd8_positive_alpha_beta_t_cell.name}
)
```

Re-run validation:

```python
try:
    curator.validate()
except ln.errors.ValidationError as error:
    print(error)
```

For the `perturbation` feature, we need to register new perturbations:

```python
ln.Record.from_values(["DMSO", "IFNG"], create=True).save()
```

### (6) Save a validated & annotated dataset

```python
artifact = curator.save_artifact(key="examples/my_curated_dataset.parquet")
artifact.describe()
```

## Common fixes

This section covers the most frequent curation issues and their solutions.
Use this as a reference when validation fails.

### Feature validation errors

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

### Value validation errors

<!-- #region -->

**Issue**: "Terms not validated in feature 'cell_type'"

```
2 terms not validated in feature 'cell_type': 'B-cell', 'CD8-pos alpha-beta T cell'
    1 synonym found: "B-cell" → "B cell"
    → curate synonyms via: curator.cat.standardize("cell_type")
    for remaining terms:
    → fix typos, remove non-existent values, or create objects via:

  objects = bionty.CellType.from_values(['CD8-pos alpha-beta T cell'], field='name').save()
```

**Solutions**:

```python
# Solution 1: Use automatic standardization if given hint (handles synonyms)
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

# Solution 4: Create new legitimate objects
objects = bt.CellType.from_values(["my_new_cell_type"]).save()
```

<!-- #endregion -->

### Data type errors

<!-- #region -->

**Issue**: "Expected categorical data, got object"

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

### Organism-specific ontology errors

<!-- #region -->

**Issue**: "Terms not validated" for organism-specific ontologies like developmental stages

```
2 terms not validated in feature 'developmental_stage_ontology_id': 'MmusDv:0000142', 'MmusDv:0000022'
```

**Solution**: Specify organism-specific source in feature definition using `cat_filters`:

```python
# When defining the schema, specify the organism-specific source
mouse_source = bt.Source.filter(
    entity="bionty.DevelopmentalStage",
    organism="mouse"
).one()

schema = ln.Schema(
    features=[
        ln.Feature(
            name="developmental_stage_ontology_id",
            dtype=bt.DevelopmentalStage.ontology_id,
            cat_filters={"source": mouse_source}  # Specify organism-specific source
        )
    ],
    ...
)
```

This pattern applies to any ontology where the same registry serves multiple organisms (e.g., `DevelopmentalStage`, `Phenotype`, ...).

<!-- #endregion -->

## External data validation

Since not all metadata is always stored within the dataset itself, it is also possible to validate external metadata.
For instance, you might want to validate a separate metadata dictionary or JSON file against a schema before attaching it to your data.

```{eval-rst}
.. literalinclude:: scripts/curate_dataframe_external_features.py
   :language: python
   :caption: curate_dataframe_external_features.py
```

```python
!python scripts/curate_dataframe_external_features.py
```

## Union dtypes

Some metadata columns might validate against several registries.
This script demonstrates how to configure a feature that accepts values from multiple sources, such as allowing either a gene symbol or an Ensembl ID.

```{eval-rst}
.. literalinclude:: scripts/curate_dataframe_union_features.py
   :language: python
   :caption: curate_dataframe_union_features.py
```

```python
!python scripts/curate_dataframe_union_features.py
```

## AnnData

`AnnData`, like all other data structures that follow, is a composite structure that stores different arrays in different `slots`.

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

Under-the-hood, this uses the following built-in schema ({func}`~lamindb.examples.schemas.anndata_ensembl_gene_ids_and_valid_features_in_obs`):

```{eval-rst}
.. literalinclude:: scripts/define_schema_anndata_ensembl_gene_ids_and_valid_features_in_obs.py
   :language: python
```

This schema transposes the `var` DataFrame during curation, so that one validates and annotates the columns of `var.T`, i.e., `[ENSG00000153563, ENSG00000010610, ENSG00000170458]`.
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
valid_cell_types = curator.slots["obs"].cat.lookup(public=True)["cell_type_by_expert"]
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

Here, we demonstrate how to curate such metadata for AnnData:

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

MuData objects contain multiple modalities, each with its own structure.
The following script shows how to define and validate schemas across different modalities (e.g., RNA and ATAC) simultaneously.

```{eval-rst}
.. literalinclude:: scripts/curate_mudata.py
   :language: python
   :caption: curate_mudata.py
```

```python
!python scripts/curate_mudata.py
```

## SpatialData

For SpatialData, we need to validate annotations nested deep within the `.attrs` dictionary.
This script illustrates how to target and curate these nested metadata fields.

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

TileDB-SOMA experiments store large-scale single-cell data on disk.
Here we show how to validate the `obs` and `var` dataframes of a SOMA experiment without loading the entire dataset into memory.

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
