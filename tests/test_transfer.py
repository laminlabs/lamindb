import lamindb as ln
import pytest


# this test has to be refactored and sped up a lot
def test_transfer_from_remote_to_local():
    """Test transfer from remote to local instance."""
    import bionty as bt

    bt.Gene.filter().delete()
    bt.Organism.filter().delete()
    ln.ULabel.filter().delete()
    bt.CellType.filter().delete()

    # transfer 1st artifact
    artifact = (
        ln.Artifact.using("laminlabs/lamin-dev")
        .filter(uid="livFRRpMaOgb3y8U2mK2")
        .one()
    )

    id_remote = artifact.id
    run_remote = artifact.run
    transform_remote = artifact.transform
    created_by_remote = artifact.created_by
    storage_remote = artifact.storage
    organism_remote = artifact.organisms.get(name="human")

    artifact.save()

    # check all ids are adjusted
    assert id_remote != artifact.id
    assert run_remote != artifact.run
    assert transform_remote != artifact.transform
    assert created_by_remote.handle != artifact.created_by.handle
    assert storage_remote.uid == artifact.storage.uid
    assert storage_remote.created_at != artifact.storage.created_at
    organism = artifact.organisms.get(name="human")
    assert organism.created_at != organism_remote.created_at

    # now check that this is idempotent and we can run it again
    artifact_repeat = (
        ln.Artifact.using("laminlabs/lamin-dev")
        .filter(uid="livFRRpMaOgb3y8U2mK2")
        .one()
    )
    artifact_repeat.save()

    # now prepare a new test case
    # mimic we have an existing feature with a different uid but same name
    feature = ln.Feature.filter(name="organism").one()
    feature.uid = "existing"
    feature.save()

    # transfer 2nd artifact
    artifact2 = (
        ln.Artifact.using("laminlabs/lamin-dev")
        .filter(uid="qz35YaRk09XtYAyLvjyZ")
        .one()
    )
    artifact2.save()

    # check the feature name
    bt.settings.organism = "mouse"
    assert artifact2.organisms.get(name="mouse") == bt.settings.organism
    assert artifact.features["obs"].get(name="organism").uid == "existing"

    artifact.delete(permanent=True, storage=False)
    artifact2.delete(permanent=True, storage=False)
    ln.FeatureSet.filter().delete()
    bt.Gene.filter().delete()
    bt.Organism.filter().delete()
    ln.ULabel.filter().delete()
    bt.Disease.filter().delete()
    bt.CellLine.filter().delete()
    bt.CellType.filter().delete()
    bt.Phenotype.filter().delete()
    bt.Ethnicity.filter().delete()
    bt.ExperimentalFactor.filter().delete()
    bt.DevelopmentalStage.filter().delete()
    bt.Tissue.filter().delete()
    ln.Feature.filter().delete()


def test_using_record_organism():
    """Test passing record and organism to the using_key instance."""
    import bionty as bt

    release_110 = bt.Source.filter(
        organism="mouse", entity="bionty.Gene", version="release-110"
    ).one()
    release_110_cxg = (
        bt.Source.using("laminlabs/lamin-dev")
        .filter(organism="mouse", entity="bionty.Gene", version="release-110")
        .one()
    )
    release_112_cxg = (
        bt.Source.using("laminlabs/lamin-dev")
        .filter(organism="mouse", entity="bionty.Gene", version="release-112")
        .one()
    )

    # passing the wrong source
    inspector = bt.Gene.using("laminlabs/lamin-dev").inspect(
        ["ENSMUSG00000102862", "ENSMUSG00000084826"],
        field=bt.Gene.ensembl_gene_id,
        source=release_112_cxg,
    )
    assert len(inspector.validated) == 0

    # passing the correct source
    inspector = bt.Gene.using("laminlabs/lamin-dev").inspect(
        ["ENSMUSG00000102862", "ENSMUSG00000084826"],
        field=bt.Gene.ensembl_gene_id,
        source=release_110_cxg,
    )
    assert len(inspector.validated) == 2

    # passing the correct source but from the wrong instance
    with pytest.raises(ValueError) as error:
        inspector = bt.Gene.using("laminlabs/lamin-dev").inspect(
            ["ENSMUSG00000102862", "ENSMUSG00000084826"],
            field=bt.Gene.ensembl_gene_id,
            source=release_110,
        )
    assert (
        "source must be a bionty.Source record from instance 'laminlabs/lamin-dev'"
        in str(error.value)
    )
