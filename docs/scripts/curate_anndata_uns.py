import lamindb as ln

ln.core.datasets.mini_immuno.define_features_labels()
adata = ln.core.datasets.mini_immuno.get_dataset1(otype="AnnData")
adata.uns["study_metadata"] = {"temperature": 21.6, "experiment_id": "EXP001"}

schema = ln.Schema.get(name="Study metadata schema")
curator = ln.curators.AnnDataCurator(adata, schema)
curator.validate()
