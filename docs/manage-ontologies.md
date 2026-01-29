[![Jupyter Notebook](https://img.shields.io/badge/Source%20on%20GitHub-orange)](https://github.com/laminlabs/lamindb/blob/main/docs/manage-ontologies.ipynb)

# Manage biological ontologies

This guide shows how to manage ontologies for basic biological entities based on plugin {mod}`bionty`.

If instead you're interested in

- accessing public ontologies, see {doc}`docs:public-ontologies`
- flexible bio registries for the wetlab (a LIMS), see {class}`~lamindb.Record` and {doc}`docs:sheets`

```python
# pip install lamindb
!lamin init --storage ./test-ontologies --modules bionty
```

## Import records from public ontologies

Let's first populate our {class}`~bionty.CellType` registry with the default public ontology (Cell Ontology).

```python
import lamindb as ln
import bionty as bt

# inspect the available public ontology versions
bt.Source.to_dataframe()
```

```python
# inspect which ontology version we're about to import
bt.Source.get(entity="bionty.CellType", currently_used=True)
```

```python
# populate the database with a public ontology
bt.CellType.import_source()
```

This is now your in-house cell type ontology in which you can add & modify records as you like. It's a registry just like `Artifact` or `Record`.

```python
# all public cell types are now available in LaminDB
bt.CellType.to_dataframe()
```

```python
# let's also populate the Gene registry with human and mouse genes
bt.Gene.import_source(organism="human")
bt.Gene.import_source(organism="mouse")
```

## Access records in in-house registries

Search key words:

```python
bt.CellType.search("gamma-delta T").to_dataframe().head(2)
```

Or look up with auto-complete:

```python
cell_types = bt.CellType.lookup()
hsc_record = cell_types.hematopoietic_stem_cell
hsc_record
```

Filter by fields and relationships:

```python
gdt_cell = bt.CellType.get(ontology_id="CL:0000798", created_by__handle="testuser1")
gdt_cell
```

View the ontological hierarchy:

```python
gdt_cell.view_parents()  # pass with_children=True to also view children
```

Or access the parents and children directly:

```python
gdt_cell.parents.to_dataframe()
```

```python
gdt_cell.children.to_dataframe()
```

It is also possible to recursively query parents or children, getting direct parents (children), their parents, and so forth.

```python
gdt_cell.query_parents().to_dataframe()
```

```python
gdt_cell.query_children().to_dataframe()
```

## Construct custom hierarchies of records

You can add a child of a parent record:

```python
# register a new cell type
my_celltype = bt.CellType(name="my new T-cell subtype").save()
# specify "gamma-delta T cell" as a parent
my_celltype.parents.add(gdt_cell)

# visualize hierarchy
my_celltype.view_parents(distance=3)
```

## Create new records

When accessing datasets, one often encounters bulk references to entities that might be corrupted or standardized using different standardization schemes.

Let's consider an example based on an `AnnData` object, in the `cell_type` annotations of this `AnnData` object, we find 4 references to cell types:

```python
adata = ln.examples.datasets.anndata_with_obs()
adata.obs.cell_type.value_counts()
```

We'd like to load the corresponding records in our in-house registry to annotate a dataset.

To this end, you'll typically use {class}`~lamindb.models.CanCurate.from_values`, which will both validate & retrieve records that match the values.

```python
cell_types = bt.CellType.from_values(adata.obs.cell_type)
cell_types
```

Logging informed us that 3 cell types were validated. Since we loaded these records at the same time, we could readily use them to annotate a dataset.

:::{dropdown} What happened under-the-hood?

`.from_values()` performs the following look ups:

1. If registry records match the values, load these records
2. If values match synonyms of registry records, load these records
3. If no record in the registry matches, attempt to load records from a public ontology
4. Same as 3. but based on synonyms

No records will be returned if all 4 look ups are unsuccessful.

Sometimes, it's useful to treat validated records differently from non-validated records. Here is a way:

```
original_values = ["gut", "gut2"]
inspector = bt.Tissue.inspect(original_values)
records_from_validated_values = bt.Tissue.from_values(inspector.validated)
```

:::

Alternatively, we can retrieve records based on ontology ids:

```python
adata.obs.cell_type_id.unique().tolist()
```

```python
bt.CellType.from_values(adata.obs.cell_type_id, field=bt.CellType.ontology_id)
```

## Validate & standardize

Simple validation of an iterable of values works like so:

```python
bt.CellType.validate(["fat cell", "blood forming stem cell"])
```

Because these values don't comply with the registry, they're not validated!

You can easily convert these values to validated standardized names based on synonyms like so:

```python
bt.CellType.standardize(["fat cell", "blood forming stem cell"])
```

Alternatively, you can use `.from_values()`, which will only ever return validated records and automatically standardize under-the-hood:

```python
bt.CellType.from_values(["fat cell", "blood forming stem cell"])
```

If you are now sure what to do, use `.inspect()` to get instructions:

```python
bt.CellType.inspect(["fat cell", "blood forming stem cell"]);
```

We can also add new synonyms to a record:

```python
hsc_record.add_synonym("HSC")
```

And when we encounter this synonym as a value, it will now be standardized using synonyms-lookup, and mapped on the correct registry record:

```python
bt.CellType.standardize(["HSC"])
```

A special synonym is `.abbr` (short for abbreviation), which has its own field and can be assigned via:

```python
hsc_record.set_abbr("HSC")
```

You can create a lookup object from the `.abbr` field:

```python
cell_types = bt.CellType.lookup("abbr")
cell_types.hsc
```

The same workflow works for all of `bionty`'s registries.

## Manage ontologies across organisms

Several registries are organism-aware (has a `.organism` field), for instance, {class}`~bionty.Gene`.

In this case, API calls that interact with multi-organism registries require an `organism` argument when there's ambiguity.

For instance, when validating gene symbols:

```python
bt.Gene.validate(["TCF7", "ABC1"], organism="human")
```

In contrary, working with Ensembl Gene IDs doesn't require passing `organism`, as there's no ambiguity:

```python
bt.Gene.validate(
    ["ENSG00000000419", "ENSMUSG00002076988"], field=bt.Gene.ensembl_gene_id
)
```

When working with the same organism throughout your analysis/workflow, you can omit the `organism` argument by configuring it globally:

```python
bt.settings.organism = "mouse"
bt.Gene.from_source(symbol="Ap5b1")
```

## Track ontology versions

Under-the-hood, source ontology versions are automatically tracked for each registry:

```python
bt.Source.filter(currently_used=True).to_dataframe()
```

Each record is linked to a versioned public source (if it was created from public):

```python
hepatocyte = bt.CellType.get(name="hepatocyte")
hepatocyte.source
```

## Create records from a specific ontology version

By default, new records are imported or created from the `"currently_used"` public sources which are configured during the instance initialization, e.g.:

```python
bt.Source.filter(entity="bionty.Phenotype", currently_used=True).to_dataframe()
```

Sometimes, the default source doesn't contain the ontology term you are looking for.

You can then specify to create a record from a non-default source. For instance, we can use the `ncbitaxon` ontology:

```python
source = bt.Source.get(entity="bionty.Organism", name="ncbitaxon")
source
```

```python
# validate against the NCBI Taxonomy
bt.Organism.validate(
    ["iris setosa", "iris versicolor", "iris virginica"], source=source
)
```

```python
# since we didn't seed the Organism registry with the NCBITaxon public ontology
# we need to save the records to the database
records = bt.Organism.from_values(
    ["iris setosa", "iris versicolor", "iris virginica"], source=source
).save()

# now we can query a iris organism and view its parents and children
bt.Organism.get(name="iris").view_parents(with_children=True)
```

<!-- #region -->

## Access any Ensembl genes

Genes from all Ensembl versions and organisms can be accessed, even though they are not yet present in the `bt.Source` registry.

For instance, if you want to use `rabbit` genes from Ensembl version `release-103`:

```python

# pip install pymysql
import bionty as bt

# automatically download genes for a new organism
gene_ontology = bt.base.Gene(source="ensembl", organism="rabbit", version='release-103')

# register the new source in lamindb
gene_ontology.register_source_in_lamindb()

# now you can start using this source

# import all genes from this source to your Gene registry
source = bt.Source.get(entity="bionty.Gene", name="ensembl", organism="rabbit", version="release-103")
bt.Gene.import_source(source=source)
```

<!-- #endregion -->
