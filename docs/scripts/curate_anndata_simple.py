import lamindb as ln

ln.core.datasets.mini_immuno.define_features_labels()
adata = ln.core.datasets.mini_immuno.get_dataset1(otype="AnnData")
schema = ln.schemas.anndata.ensembl_gene_ids_and_valid_features_in_obs()
artifact = ln.Artifact.from_anndata(
    adata, key="my_datasets/mini_immuno.h5ad", schema=schema
).save()
artifact.describe()
