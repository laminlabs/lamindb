import lamindb as ln

ln.examples.datasets.mini_immuno.define_features_labels()
adata = ln.examples.datasets.mini_immuno.get_dataset1(otype="AnnData")
adata.uns["study_metadata"] = adata.uns.copy()
schema = ln.Schema.get(name="anndata_nested_study_metadata")
artifact = ln.Artifact.from_anndata(
    adata, schema=schema, key="examples/mini_immuno_uns.h5ad"
)
artifact.describe()
