import lamindb as ln
import bionty as bt

# define valid labels
perturbation = ln.ULabel(name="Perturbation", is_type=True).save()
ln.ULabel(name="DMSO", type=perturbation).save()
ln.ULabel(name="IFNG", type=perturbation).save()
bt.CellType.from_source(name="B cell").save()
bt.CellType.from_source(name="T cell").save()

# define obs schema
obs_schema = ln.Schema(
    name="small_dataset1_obs_level_metadata",
    features=[
        ln.Feature(name="perturbation", dtype="cat[ULabel[Perturbation]]").save(),
        ln.Feature(name="sample_note", dtype=str).save(),
        ln.Feature(name="cell_type_by_expert", dtype=bt.CellType).save(),
        ln.Feature(name="cell_type_by_model", dtype=bt.CellType).save(),
    ],
).save()

# define var schema
var_schema = ln.Schema(
    name="scRNA_seq_var_schema",
    itype=bt.Gene.ensembl_gene_id,
    dtype=int,
).save()

# define composite schema
anndata_schema = ln.Schema(
    name="small_dataset1_anndata_schema",
    otype="AnnData",
    components={"obs": obs_schema, "var": var_schema},
).save()

# curate an AnnData
adata = ln.core.datasets.small_dataset1(otype="AnnData")
curator = ln.curators.AnnDataCurator(adata, anndata_schema)
artifact = curator.save_artifact(key="example_datasets/dataset1.h5ad")
assert artifact.schema == anndata_schema
