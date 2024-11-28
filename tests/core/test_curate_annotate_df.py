import anndata as ad
import bionty as bt
import lamindb as ln
import pandas as pd
from lamindb.core._data import _describe_postgres


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

    # define the data in the dataset
    # it's a mix of numerical measurements and observation-level metadata
    dataset_dict = {
        "CD8A": [1, 2, 3],
        "CD4": [3, 4, 5],
        "CD14": [5, 6, 7],
        "medium": ["DMSO", "IFNG", "DMSO"],
        "sample_note": ["was ok", "looks naah", "pretty! 🤩"],
        "cell_type_by_expert": ["B cell", "T cell", "T cell"],
        "cell_type_by_model": ["B cell", "T cell", "T cell"],
    }
    # define the dataset-level metadata
    metadata = {
        "temperature": 21.6,
        "study": "Candidate marker study 1",
        "date_of_experiment": "2024-12-01",
        "experiment_note": "We had a great time performing this experiment and the results look compelling.",
    }
    # the dataset as DataFrame
    dataset_df = pd.DataFrame(dataset_dict, index=["sample1", "sample2", "sample3"])
    dataset_ad = ad.AnnData(dataset_df.iloc[:, :3], obs=dataset_df.iloc[:, 3:])
    # curate dataset
    curator = ln.Curator.from_anndata(
        dataset_ad,
        var_index=bt.Gene.symbol,
        categoricals={
            "medium": ln.ULabel.name,
            "cell_type_by_expert": bt.CellType.name,
            "cell_type_by_model": bt.CellType.name,
        },
        organism="human",
    )
    artifact = curator.save_artifact(key="example_datasets/dataset1.h5ad")
    # annotate with dataset-level features
    artifact.features.add_values(metadata)

    # expected output has italicized elements that can't be tested
    # hence testing is restricted to section content, not headings
    description = _describe_postgres(artifact, print_types=True)
    # > Artifact(uid='tQAwzih2n44VQRjO0000', is_latest=True, key='example_datasets/dataset1.h5ad', suffix='.h5ad', type='dataset', size=23560, hash='voB-uoihaivmNskhV7osPQ', n_observations=3, _hash_type='md5', _accessor='AnnData', visibility=1, _key_is_virtual=True, created_at=2024-11-28 17:05:30 UTC)
    # >   Provenance
    # >     .storage: Storage = '/Users/falexwolf/repos/laminhub/rest-hub/sub/lamindb/default_storage_unit_core'
    # >     .created_by: User = 'falexwolf'
    # >   Labels
    # >     .cell_types: CellType = 'B cell', 'T cell'
    # >     .ulabels: ULabel = 'DMSO', 'IFNG', 'Candidate marker study 1'
    # >   Feature sets
    # >     'var' = 'CD8A', 'CD4', 'CD14'
    # >     'obs' = 'medium', 'sample_note', 'cell_type_by_expert', 'cell_type_by_model'
    # >   Feature values -- internal
    # >     'cell_type_by_expert': cat[bionty.CellType] = B cell, T cell
    # >     'cell_type_by_model': cat[bionty.CellType] = B cell, T cell
    # >     'medium': cat[ULabel] = DMSO, IFNG
    # >   Feature values -- external
    # >     'study': cat[ULabel] = Candidate marker study 1
    # >     'date_of_experiment': date = 2024-12-01
    # >     'experiment_note': str = We had a great time performing this experiment and the results look compelling.
    # >     'temperature': float = 21.6

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
