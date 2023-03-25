# The story of data objects in story & memory

We have come to love the pydata family of `DataFrame`, `AnnData`, `pytorch.DataLoader`, `zarr.Array`, `pyarrow.Table`, `xarray.Dataset`, and others for accessing lower-level data objects.

But we couldn’t find an object for accessing how data objects are linked to context.
So, we made `File` to help with modeling and understanding data objects in relation to their context.

Context can be other data objects, data transformations, ML models, people & pipelines who performed transformations, and all aspects of data lineage.
Context can also be hypotheses and any entity of the domain in which data is generated and modeled.

Depending on how `File`s are linked to context, they give rise to features of data lakes, warehouses and knowledge graphs.

We focused on linking `File` to biological concepts: entities, their types, records, transformations, and relations.
You'll learn about them further down the guide.

```{Note}

Learn more: {class}`lamindb.File`.
```
