import lamindb as ln


def test_save_ulabel_to_another_db_via_model_save():
    assert ln.setup.settings.instance.name == "testdb2"

    using = f"{ln.setup.settings.user.handle}/testdb1"
    db2 = ln.DB(using)
    name = "save-using-single"
    remote_labels = db2.ULabel.filter(name=name)
    if remote_labels.exists():
        remote_labels.delete(permanent=True)

    try:
        ulabel = ln.ULabel(name=name).save(using=using)

        assert ulabel._state.db == using
        assert db2.ULabel.get(name=name)
        assert ln.ULabel.filter(name=name).count() == 0
    finally:
        db2.ULabel.filter(name=name).delete(permanent=True)


def test_save_ulabels_to_another_db_via_ln_save():
    assert ln.setup.settings.instance.name == "testdb2"

    using = f"{ln.setup.settings.user.handle}/testdb1"
    db2 = ln.DB(using)
    names = [f"save-using-bulk-{i}" for i in range(3)]
    remote_labels = db2.ULabel.filter(name__in=names)
    if remote_labels.exists():
        remote_labels.delete(permanent=True)

    try:
        ulabels = [ln.ULabel(name=name) for name in names]
        ln.save(ulabels, using=using)

        remote_matches = db2.ULabel.filter(name__in=names)
        assert remote_matches.count() == len(names)
        assert ln.ULabel.filter(name__in=names).count() == 0
    finally:
        db2.ULabel.filter(name__in=names).delete(permanent=True)


def test_save_record_with_relation_to_another_db_via_model_save():
    assert ln.setup.settings.instance.name == "testdb2"

    using = f"{ln.setup.settings.user.handle}/testdb1"
    db2 = ln.DB(using)
    type_name = "save-using-relation-type"
    rec_name = "save-using-relation-child"

    for qs in (db2.Record.filter(name=rec_name), db2.Record.filter(name=type_name)):
        if qs.exists():
            qs.delete(permanent=True)

    try:
        rec_type = ln.Record(name=type_name, is_type=True).save(using=using)
        child = ln.Record(name=rec_name, type=rec_type).save(using=using)

        assert child._state.db == using
        assert db2.Record.get(name=rec_name)
        assert ln.Record.filter(name=rec_name).count() == 0
    finally:
        for qs in (db2.Record.filter(name=rec_name), db2.Record.filter(name=type_name)):
            qs.delete(permanent=True)
