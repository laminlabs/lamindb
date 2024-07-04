import lamindb as ln


# this test has to be refactored and sped up a lot
def test_transfer_from_remote_to_local():
    """Test transfer from remote to local instance."""
    import bionty as bt

    bt.Gene.filter().delete()
    bt.Organism.filter().delete()
    ln.ULabel.filter().delete()

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

    artifact.save(parents=False)

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
    artifact_repeat.save(parents=False)

    # now prepare a new test case
    # mimic we have an existing feature with a different uid but same name
    feature = ln.Feature.filter(name="organism").one()
    feature.uid = "existing"
    feature.save()

    # transfer 2nd artifact
    bt.settings.auto_save_parents = False
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
