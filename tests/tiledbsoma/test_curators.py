import shutil

import bionty as bt
import lamindb as ln
import pytest
import tiledbsoma
import tiledbsoma.io


def test_tiledbsoma_curator(clean_soma_files):
    """Test TiledbSomaExperimentCurator with schema."""
    obs_schema = ln.Schema(
        features=[
            ln.Feature(name="cell_type_by_expert", dtype=bt.CellType).save(),
            ln.Feature(name="cell_type_by_model", dtype=bt.CellType).save(),
        ],
    ).save()

    var_schema = ln.Schema(
        features=[
            ln.Feature(name="var_id", dtype=bt.Gene.ensembl_gene_id).save(),
        ],
        coerce=True,
    ).save()

    soma_schema = ln.Schema(
        otype="tiledbsoma",
        slots={
            "obs": obs_schema,
            "ms:RNA": var_schema,
        },
    ).save()

    # Convert AnnData to SOMA format
    adata = ln.examples.datasets.mini_immuno.get_dataset1(otype="AnnData")
    tiledbsoma.io.from_anndata(
        "small_dataset.tiledbsoma", adata, measurement_name="RNA"
    )

    # Test with invalid dataset
    with pytest.raises(ln.errors.InvalidArgument) as e:
        ln.curators.TiledbsomaExperimentCurator(adata, soma_schema)
    assert "dataset must be SOMAExperiment-like." in str(e.value)

    # Test with invalid schema
    with tiledbsoma.Experiment.open("small_dataset.tiledbsoma") as experiment:
        with pytest.raises(ln.errors.InvalidArgument) as e:
            ln.curators.TiledbsomaExperimentCurator(experiment, schema=var_schema)
        assert "Schema otype must be 'tiledbsoma'." in str(e.value)

    with tiledbsoma.Experiment.open("small_dataset.tiledbsoma") as experiment:
        curator = ln.curators.TiledbsomaExperimentCurator(experiment, soma_schema)

        assert "obs" in curator.slots
        assert "ms:RNA" in curator.slots

        curator.validate()

        artifact = curator.save_artifact(
            key="examples/soma_experiment.tiledbsoma",
            description="SOMA experiment with schema validation",
        )

        assert artifact.schema == soma_schema
        assert "obs" in artifact.features.slots
        assert "ms:RNA" in artifact.features.slots

        # Check feature values are properly annotated
        assert set(artifact.features.get_values()["cell_type_by_expert"]) == {
            "CD8-positive, alpha-beta T cell",
            "B cell",
        }
        assert set(artifact.features.get_values()["cell_type_by_model"]) == {
            "T cell",
            "B cell",
        }

    # Altered data (gene typo)
    adata_typo = ln.examples.datasets.mini_immuno.get_dataset1(
        otype="AnnData", with_gene_typo=True
    )
    typo_soma_path = "./mini_immuno_dataset1_typo.tiledbsoma"
    tiledbsoma.io.from_anndata(typo_soma_path, adata_typo, measurement_name="RNA")
    with tiledbsoma.Experiment.open(typo_soma_path) as experiment_typo:
        curator_typo = ln.curators.TiledbsomaExperimentCurator(
            experiment_typo, soma_schema
        )

        # Validation should fail due to typo
        with pytest.raises(ln.errors.ValidationError) as error:
            curator_typo.validate()
        assert "GeneTypo" in str(error.value)

    # Clean up
    shutil.rmtree(typo_soma_path)
    artifact.delete(permanent=True)
    soma_schema.delete(permanent=True)
    var_schema.delete(permanent=True)
    obs_schema.delete(permanent=True)
