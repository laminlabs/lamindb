import anndata as ad
import bionty as bt
import lamindb as ln
import pandas as pd
from lamindb.core._data import describe


def test_curate_annotate_df():
    ## Define the schema of the dataset & its metadata

    # observation-level metadata
    ln.Feature(name="medium", dtype="cat[ULabel]").save()
    ln.Feature(name="sample_note", dtype="str").save()
    ln.Feature(name="cell_type_by_expert", dtype="cat[bionty.CellType]").save()
    ln.Feature(name="cell_type_by_model", dtype="cat[bionty.CellType]").save()
    # dataset-level metadata
    ln.Feature(name="temperature", dtype="float").save()
    ln.Feature(name="study", dtype="cat[ULabel]").save()
    ln.Feature(name="date_of_experiment", dtype="date").save()
    ln.Feature(name="experiment_note", dtype="str").save()

    ## Register permissible values for categoricals

    ln.save(ln.ULabel.from_values(["DMSO", "IFNG"], create=True))
    ln.save(
        ln.ULabel.from_values(
            ["Candidate marker study 1", "Candidate marker study 2"], create=True
        )
    )
    ln.save(bt.CellType.from_values(["B cell", "T cell"], create=True))

    ## Ingest a dataset

    dataset = pd.DataFrame(
        {
            "CD8A": [1, 2, 3],
            "CD4": [3, 4, 5],
            "CD14": [5, 6, 7],
            "medium": ["DMSO", "IFNG", "DMSO"],
            "sample_note": ["was ok", "looks naah", "pretty! 🤩"],
            "cell_type_by_expert": ["B cell", "T cell", "T cell"],
            "cell_type_by_model": ["B cell", "T cell", "T cell"],
        },
        index=["sample1", "sample2", "sample3"],
    )
    metadata = {
        "temperature": 21.6,
        "study": "Candidate marker study 1",
        "date_of_experiment": "2024-12-01",
        "experiment_note": "We had a great time performing this experiment and the results look compelling.",
    }
    adata1 = ad.AnnData(dataset.iloc[:, :3], obs=dataset.iloc[:, 3:])
    # curate dataset
    curator = ln.Curator.from_anndata(
        adata1,
        var_index=bt.Gene.symbol,
        categoricals={
            "medium": ln.ULabel.name,
            "cell_type_by_expert": bt.CellType.name,
            "cell_type_by_model": bt.CellType.name,
        },
        organism="human",
    )
    artifact1 = curator.save_artifact(key="example_datasets/dataset1.h5ad")
    # annotate with dataset-level features
    artifact1.features.add_values(metadata)

    description = describe(artifact1, print_types=True)
    print(description)
    labels = """.cell_types: bionty.CellType = 'B cell', 'T cell'
    .ulabels: ULabel = 'DMSO', 'IFNG', 'Candidate marker study 1'"""
    assert labels in description

    internal_features = """'cell_type_by_expert': cat[bionty.CellType] = B cell, T cell
    'cell_type_by_model': cat[bionty.CellType] = B cell, T cell
    'medium': cat[ULabel] = DMSO, IFNG"""
    assert internal_features in description

    external_features = """'study': cat[ULabel] = Candidate marker study 1
    'date_of_experiment': date = 2024-12-01
    'experiment_note': str = We had a great time performing this experiment and the results look compelling.
    'temperature': float = 21.6"""
    assert external_features in description
