import lamindb as ln

# define features and labels
ln.core.datasets.mini_immuno.define_features_labels()

# define a schema
schema = ln.schemas.anndata.ensembl_gene_ids_and_valid_features_in_obs()

# curate an AnnData
adata = ln.core.datasets.mini_immuno.get_dataset1(otype="AnnData")
artifact = ln.Artifact.from_anndata(adata, schema=schema).save()
assert artifact.schema == schema
assert artifact.features.slots["var.T"].members.df()["ensembl_gene_id"].tolist() == [
    "ENSG00000153563",
    "ENSG00000010610",
    "ENSG00000170458",
]
