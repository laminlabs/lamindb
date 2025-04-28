import lamindb as ln
import bionty as bt


# define features and labels
ln.core.datasets.mini_immuno.define_features_labels()

# define a schema
obs_schema = ln.Schema(itype=ln.Feature).save()
varT_schema = ln.Schema(
    itype=bt.Gene.ensembl_gene_id
).save()  # map columns on Gene registry
schema = ln.Schema(
    otype="AnnData",
    components={"obs": obs_schema, "var.T": varT_schema},
).save()

# curate an AnnData
adata = ln.core.datasets.mini_immuno.get_dataset1(otype="AnnData")
artifact = ln.Artifact.from_anndata(adata, schema=schema).save()
assert artifact.schema == schema
