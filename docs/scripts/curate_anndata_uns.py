import lamindb as ln

ln.examples.datasets.mini_immuno.define_features_labels()
adata = ln.examples.datasets.mini_immuno.get_dataset1(otype="AnnData")
schema = ln.Schema.get(name="Study metadata schema")
artifact = ln.Artifact.from_anndata(
    adata, schema=schema, key="examples/mini_immuno_uns.h5ad"
)
artifact.describe()
