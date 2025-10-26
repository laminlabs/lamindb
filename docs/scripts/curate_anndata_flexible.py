import lamindb as ln

ln.examples.datasets.mini_immuno.define_features_labels()
adata = ln.examples.datasets.mini_immuno.get_dataset1(otype="AnnData")
artifact = ln.Artifact.from_anndata(
    adata,
    key="examples/mini_immuno.h5ad",
    schema="ensembl_gene_ids_and_valid_features_in_obs",
).save()
artifact.describe()
