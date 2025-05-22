import lamindb as ln
import bionty as bt

# observation-level metadata
ln.Feature(name="perturbation", dtype="cat[ULabel]").save()
ln.Feature(name="sample_note", dtype="str").save()
ln.Feature(name="cell_type_by_expert", dtype="cat[bionty.CellType]").save()
ln.Feature(name="cell_type_by_model", dtype="cat[bionty.CellType]").save()
# dataset-level metadata
ln.Feature(name="temperature", dtype="float").save()
ln.Feature(name="experiment", dtype="cat[ULabel]").save()
ln.Feature(name="date_of_study", dtype="date").save()
ln.Feature(name="study_note", dtype="str").save()
ln.Feature(name="study_metadata", dtype=dict).save()

## Permissible values for categoricals
ln.ULabel.from_values(["DMSO", "IFNG"], create=True).save()
ln.ULabel.from_values(["Experiment 1", "Experiment 2"], create=True).save()
bt.CellType.from_values(["B cell", "T cell"], create=True).save()

schema = ln.examples.schemas.anndata_ensembl_gene_ids_and_valid_features_in_obs()

## Ingest dataset1
adata = ln.core.datasets.mini_immuno.get_dataset1(otype="AnnData")
artifact = ln.Artifact.from_anndata(
    adata,
    key="examples/dataset1.h5ad",
    schema=schema,
).save()
adhoc = {"study_metadata": {"detail1": "123", "detail2": 1}}
dataset_metadata = adata.uns
dataset_metadata.update(adhoc)
artifact.features.add_values(dataset_metadata)  # type: ignore

# Ingest dataset2
adata2 = ln.core.datasets.mini_immuno.get_dataset2(otype="AnnData")
artifact2 = ln.Artifact.from_anndata(
    adata2,
    key="examples/dataset2.h5ad",
    schema=schema,
).save()
adhoc2 = {"study_metadata": {"detail1": "456", "detail2": 2}}
dataset_metadata2 = adata2.uns
dataset_metadata2.update(adhoc2)
artifact2.features.add_values(dataset_metadata2)  # type: ignore
