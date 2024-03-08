import lamindb as ln


def test_delete_record():
    names = ["label1", "label2", "label3"]
    labels = [ln.ULabel(name=name) for name in names]
    ln.save(labels)
    ln.ULabel.filter(name__in=names).delete()
    assert ln.ULabel.filter(name__in=names).count() == 0
