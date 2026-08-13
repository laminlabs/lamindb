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

def test_save_record_with_multivalued_relation_to_another_db_via_model_save():
    assert ln.setup.settings.instance.name == "testdb2"

    using = f"{ln.setup.settings.user.handle}/testdb1"
    db2 = ln.DB(using)

    type_name = "mv-intervention-type"
    gene_type_name = "mv-gene-type"
    gene_feat_name = "mv-genes"
    gene_names = [f"mv-gene-{i}" for i in range(5)]
    rec_name = "mv-child"

    def _clean():
        for qs in (
            db2.Record.filter(name=rec_name),
            db2.Record.filter(name__in=gene_names),
            db2.Record.filter(name=type_name),
            db2.Record.filter(name=gene_type_name),
            db2.Feature.filter(name=gene_feat_name),
        ):
            if qs.exists():
                qs.delete(permanent=True)

    _clean()
    try:
        intervention_t = ln.Record(name=type_name, is_type=True).save(using=using)
        gene_t = ln.Record(name=gene_type_name, is_type=True).save(using=using)
        gene_feat = ln.Feature(name=gene_feat_name, dtype=gene_t).save(using=using)
        gene_pool = [ln.Record(name=g, type=gene_t).save(using=using) for g in gene_names]

        child = ln.Record(
            name=rec_name,
            type=intervention_t,
            features={gene_feat: gene_pool[:3]},
        ).save(using=using)

        assert child._state.db == using
        assert db2.Record.get(name=rec_name)
        assert ln.Record.filter(name=rec_name).count() == 0
    finally:
        _clean()
