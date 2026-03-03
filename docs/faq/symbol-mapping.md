---
execute_via: python
---

# Why should I not index datasets with gene symbols?

Gene symbols are widely used for readability, particularly for visualization. However, indexing datasets with gene symbols presents challenges:

- A single gene may have multiple symbols or aliases.
- Gene symbols change over time (e.g., _BRCA2_ was once _FACD_) without version tracking.
- The same symbol can represent different genes across species.
- Symbols may be misinterpreted by software (e.g., _SEPT9_ as "September 9" in Excel).
- Formatting inconsistencies exist (e.g., case sensitivity, special characters).

Using unique identifiers like ENSEMBL gene IDs addresses these issues by providing:

- A direct, stable mapping to genomic coordinates.
- Consistency across databases.
- Species-specific prefixes to prevent cross-species confusion.
- Unique, permanent identifiers with standardized formatting.

Storing ENSEMBL gene IDs alongside gene symbols offers readability for visualization while maintaining robust data integrity. During curation, validating against ENSEMBL gene IDs ensures accurate mapping.

If only symbols are available for a dataset, you can map them to ENSEMBL IDs using {meth}`~bionty.Gene.standardize`.

```python
# !pip install 'lamindb[bionty]'
!lamin init --storage test-symbols --modules bionty
```

```python
import lamindb as ln
import bionty as bt
import numpy as np
import pandas as pd
import anndata as ad

# create example AnnData object with gene symbols
rng = np.random.default_rng(42)
X = rng.integers(0, 100, size=(5, 10))
var = pd.DataFrame(
    index=pd.Index(
        [
            "BRCA1",
            "TP53",
            "EGFR",
            "KRAS",
            "PTEN",
            "MYC",
            "VEGFA",
            "IL6",
            "TNF",
            "GAPDH",
        ],
        name="symbol",
    )
)
adata = ad.AnnData(X=X, var=var)
adata.var
```

```python
# map Gene symbols to ENSEMBL IDs
gene_mapper = bt.Gene.standardize(
    adata.var.index,
    field=bt.Gene.symbol,
    return_field=bt.Gene.ensembl_gene_id,
    return_mapper=True,
    organism="human",
)
adata.var["ensembl_id"] = adata.var.index.map(
    lambda gene_id: gene_mapper.get(gene_id, gene_id)
)
adata.var
```

```python
standardized_genes = bt.Gene.from_values(
    [
        "ENSG00000141510",
        "ENSG00000133703",
        "ENSG00000111640",
        "ENSG00000171862",
        "ENSG00000204490",
        "ENSG00000112715",
        "ENSG00000146648",
        "ENSG00000136997",
        "ENSG00000012048",
        "ENSG00000136244",
    ],
    field=bt.Gene.ensembl_gene_id,
    organism="human",
)
ln.save(standardized_genes)
```

This allows for validating the the `ensembl_id` against the `Gene` registry using the `bt.Gene.ensembl_gene_id` field.

```python
bt.Gene.validate(adata.var["ensembl_id"], field=bt.Gene.ensembl_gene_id)
```

```{note}
Gene symbols do not map one-to-one with ENSEMBL IDs. A single gene symbol may correspond to multiple ENSEMBL IDs due to:

1. **Gene Paralogs**: Similar symbols can be shared among paralogous genes within the same species, resulting in one symbol linking to multiple ENSEMBL IDs.
2. **Pseudogenes**: Some symbols represent both functional genes and their non-functional pseudogenes, each with distinct ENSEMBL IDs.
3. **Transcript Variants**: One symbol may map to multiple ENSEMBL transcript IDs, each representing different isoforms or splice variants.

{meth}`~bionty.Gene.standardize` retrieves the first match in cases of multiple hits, which is generally sufficient but not perfectly accurate.
```

```python
!lamin delete --force test-symbols
```
