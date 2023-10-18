import lamindb as ln


def test_transfer():
    import lnschema_bionty as lb

    lb.Gene.filter().all().delete()
    lb.Species.filter().all().delete()
    ln.ULabel.filter().all().delete()

    # insert human as species id=2
    lb.settings.species = "mouse"
    lb.settings.species = "human"

    # transfer 1st file
    file = ln.File.filter(
        using="sunnyosun/cellxgene-census",
        description__icontains="tabula sapiens - lung",
    ).one()
    ulabel_id_remote = file.ulabels.get(name="Tabula Sapiens").id
    file.save()

    assert file.species.get(name="human").id == lb.settings.species.id
    assert file.ulabels.get(name="Tabula Sapiens").id != ulabel_id_remote
    ulabel = file.ulabels.get(name="Tabula Sapiens")
    # mimic we have an existing ulabel with a different uid but same name
    ulabel.uid = "existing"
    ulabel.save()

    # transfer 2nd file
    file = ln.File.filter(
        using="sunnyosun/cellxgene-census",
        description__icontains="tabula sapiens - liver",
    ).one()
    file.save()

    assert file.species.get(name="human").id == lb.settings.species.id
    assert file.ulabels.get(name="Tabula Sapiens").id != ulabel_id_remote

    lb.Gene.filter().delete()
    lb.Species.filter().delete()
    ln.ULabel.filter().delete()
    lb.Disease.filter().delete()
    lb.CellLine.filter().delete()
    ln.CellType.filter().delete()
