import lamindb as ln


# this test has to be refactored and sped up a lot
def test_transfer_from_remote_to_local():
    """Test transfer from remote to local instance."""
    import bionty as bt

    bt.Gene.filter().delete()
    bt.Organism.filter().delete()
    ln.ULabel.filter().delete()

    bt.settings.organism = "human"

    # transfer 1st artifact
    artifact = (
        ln.Artifact.using("laminlabs/cellxgene")
        .filter(
            description__icontains="tabula sapiens",
        )
        .first()
    )

    id_remote = artifact.id
    run_remote = artifact.run
    transform_remote = artifact.transform
    created_by_remote = artifact.created_by
    storage_remote = artifact.storage
    organism_remote = artifact.organism.get(name="human")

    artifact.save(parents=False)

    # check all ids are adjusted
    assert artifact.organism.get(name="human") == bt.settings.organism
    assert id_remote != artifact.id
    assert run_remote != artifact.run
    assert transform_remote != artifact.transform
    assert created_by_remote.handle != artifact.created_by.handle
    assert storage_remote.uid == artifact.storage.uid
    assert storage_remote.created_at != artifact.storage.created_at
    organism = artifact.organism.get(name="human")
    assert organism != organism_remote

    # now check that this is idempotent and we can run it again
    artifact_repeat = (
        ln.Artifact.using("laminlabs/cellxgene")
        .filter(
            description__icontains="tabula sapiens",
        )
        .first()
    )
    artifact_repeat.save(parents=False)

    # now prepare a new test case
    # mimic we have an existing feature with a different uid but same name
    feature = ln.Feature.filter(name="organism").one()
    feature.uid = "existing"
    feature.save()

    # transfer 2nd artifact
    artifact2 = (
        ln.Artifact.using("laminlabs/cellxgene")
        .filter(
            description__icontains="tabula sapiens",
        )
        .last()
    )
    artifact2.save(parents=False)

    assert artifact2.organism.get(name="human") == bt.settings.organism
    assert artifact.features["obs"].get(name="organism").uid == "existing"

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
    ln.FeatureSet.filter().delete()
    ln.Run.filter().delete()
    # ln.Transform.filter().delete()
    ln.Artifact.filter().delete()
