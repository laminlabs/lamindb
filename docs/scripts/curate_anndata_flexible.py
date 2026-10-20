import lamindb as ln

ln.examples.datasets.mini_immuno.define_features_labels()
adata = ln.examples.datasets.mini_immuno.get_dataset1(otype="AnnData")
schema = ln.Schema.get(name="anndata_ensembl_gene_ids_and_obs_feature_type")
artifact = ln.Artifact.from_anndata(
    adata,
    key="examples/mini_immuno.h5ad",
    schema=schema,
).save()
artifact.describe()
