import lamindb as ln


def test_transfer():
    import lnschema_bionty as lb

    lb.Gene.filter().delete()
    lb.Organism.filter().delete()
    ln.ULabel.filter().delete()

    lb.settings.organism = "human"

    # transfer 1st file
    file = (
        ln.File.using("laminlabs/cellxgene-census")
        .filter(
            description__icontains="tabula sapiens - lung",
        )
        .one()
    )

    id_remote = file.id
    run_remote = file.run
    transform_remote = file.transform
    created_by_remote = file.created_by
    storage_remote = file.storage
    ulabel_remote = file.ulabels.get(name="Tabula Sapiens")

    file.save()

    # check all ids are adjusted
    assert file.organism.get(name="human") == lb.settings.organism
    assert id_remote != file.id
    assert run_remote != file.run
    assert transform_remote != file.transform
    assert created_by_remote.handle != file.created_by.handle
    assert storage_remote.uid == file.storage.uid
    assert storage_remote.created_at != file.storage.created_at
    ulabel = file.ulabels.get(name="Tabula Sapiens")
    assert ulabel != ulabel_remote
    # mimic we have an existing ulabel with a different uid but same name
    ulabel.uid = "existing"
    ulabel.save()

    # transfer 2nd file
    file2 = (
        ln.File.using("laminlabs/cellxgene-census")
        .filter(
            description__icontains="tabula sapiens - liver",
        )
        .one()
    )
    file2.save()

    assert file2.organism.get(name="human") == lb.settings.organism
    assert file2.ulabels.get(name="Tabula Sapiens").uid == "existing"

    lb.Gene.filter().delete()
    lb.Organism.filter().delete()
    ln.ULabel.filter().delete()
    lb.Disease.filter().delete()
    lb.CellLine.filter().delete()
    lb.CellType.filter().delete()
    ln.Run.filter().delete()
    ln.Transform.filter().delete()
    ln.File.filter().delete()
