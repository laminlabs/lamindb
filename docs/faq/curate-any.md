---
execute_via: python
---

# How do I validate & annotate arbitrary data structures?

This guide walks through the low-level API that lets you validate iterables.

You can then use the records create inferred during validation to annotate a dataset.

:::{dropdown} How do I validate based on a public ontology?

LaminDB makes it easy to validate categorical variables based on registries that inherit from {class}`~lamindb.models.CanCurate`.

{class}`~lamindb.models.CanCurate` methods validate against the registries in your LaminDB instance.
In {doc}`/manage-ontologies`, you'll see how to extend standard validation to validation against _public references_ using a `PubliOntology` object, e.g., via `public_genes = bt.Gene.public()`.
By default, {meth}`~lamindb.models.CanCurate.from_values` considers a match in a public reference a validated value for any {mod}`bionty` entity.

:::

```python
# pip install 'lamindb[zarr]'
!lamin init --storage ./test-curate-any --modules bionty
```

Define a test dataset.

```python
import lamindb as ln
import bionty as bt
import zarr
import numpy as np

data = zarr.open_group(store="data.zarr", mode="a")

data.create_dataset(name="temperature", shape=(3,), dtype="float32")
data.create_dataset(name="knockout_gene", shape=(3,), dtype=str)
data.create_dataset(name="disease", shape=(3,), dtype=str)

data["knockout_gene"][:] = np.array(
    ["ENSG00000139618", "ENSG00000141510", "ENSG00000133703"]
)
data["disease"][:] = np.random.default_rng().choice(
    ["MONDO:0004975", "MONDO:0004980"], 3
)
```

## Validate and standardize vectors

Read the `disease` array from the zarr group into memory.

```python
disease = data["disease"][:]
```

{meth}`~lamindb.models.CanCurate.validate` validates vectore-like values against reference values in a registry.
It returns a boolean vector indicating where a value has an exact match in the reference values.

```python
bt.Disease.validate(disease, field=bt.Disease.ontology_id)
```

When validation fails, you can call {meth}`~lamindb.models.CanCurate.inspect` to figure out what to do.

{meth}`~lamindb.models.CanCurate.inspect` applies the same definition of validation as {meth}`~lamindb.models.CanCurate.validate`, but returns a rich return value {class}`~lamindb.models.InspectResult`. Most importantly, it logs recommended curation steps that would render the data validated.

Note: you can use {meth}`~lamindb.models.CanCurate.standardize` to standardize synonyms.

```python
bt.Disease.inspect(disease, field=bt.Disease.ontology_id)
```

Bulk creating records using {meth}`~lamindb.models.CanCurate.from_values` only returns validated records.

```python
diseases = bt.Disease.from_values(disease, field=bt.Disease.ontology_id).save()
```

Repeat the process for more labels:

```python
experiments = ln.Record.from_values(
    ["Experiment A", "Experiment B"],
    field=ln.Record.name,
    create=True,  # create non-validated labels
).save()
genes = bt.Gene.from_values(
    data["knockout_gene"][:], field=bt.Gene.ensembl_gene_id
).save()
```

## Annotate the dataset

Register the dataset as an artifact:

```python
artifact = ln.Artifact("data.zarr", key="my_dataset.zarr").save()
```

Annotate with features:

```python
ln.Feature(name="experiment", dtype=ln.Record).save()
ln.Feature(name="disease", dtype=bt.Disease.ontology_id).save()
ln.Feature(name="knockout_gene", dtype=bt.Gene.ensembl_gene_id).save()
artifact.features.add_values(
    {"experiment": experiments, "knockout_gene": genes, "disease": diseases}
)
artifact.describe()
```

```python
# clean up test instance
!rm -r data.zarr
!rm -r ./test-curate-any
!lamin delete --force test-curate-any
```
