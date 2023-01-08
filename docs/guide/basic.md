# Track data

If you are dumping data in local or cloud storage but not able to find them easily, try LaminDB as it comes with tracking data lineage and edit history.

It takes 5min to [set up your LaminDB instance](https://lamin.ai/docs/db/guide/init), start [ingesting](https://lamin.ai/docs/db/guide/ingest) data and tracking them all!

## Data objects in storage & memory: `DObject`

{class}`~lamindb.DObject` allows accessing atomic data in object storage and loading or streaming them into memory.

## On-disk and in-memory representations of data

Data objects often have **canonical on-disk and in-memory representations**. If choices among these representations are made, a one-to-one mapping can be achieved, which means that any given dobject has a default in-memory and on-disk representation.

LaminDB offers meaningful default choices. For instance,

- It defaults to pandas DataFrames for in-memory representation of tables and allows you to configure loading tables into polars DataFrames.
- It defaults to the .parquet format for tables, but allows you to configure .csv or .ipc.
- Some datasets do not have a canonical in-memory representation, for instance, .fastq, .vcf, or files describing QC of datasets.

```{tip}

Examples for storage ⟷ memory correspondence:

- Table: .csv, .tsv, .parquet, .ipc (.feather) ⟷ pandas.DataFrame, polars.DataFrame
- Annotated matrix: .h5ad, .h5mu, .zarrad ⟷ anndata.AnnData, mudata.MuData
- Image: .jpg, .png ⟷ np.ndarray, or a dedicated imaging in-memory container
- Tensor: zarr directory, TileDB store ⟷ zarr loader, TileDB loader
- Fastq: .fastq ⟷ /
- VCF: .vcf ⟷ /
- QC: .html ⟷ /

```

```{toctree}
:maxdepth: 1
:hidden:

ingest
select
add-delete
view
```
