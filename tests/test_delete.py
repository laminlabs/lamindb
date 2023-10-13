import lamindb as ln


def test_delete():
    names = ["label1", "label2", "label3"]
    labels = [ln.ULabel(name=name) for name in names]
    ln.save(labels)
    ln.delete(labels[0])
    ln.delete([labels[1], labels[2]])
    assert ln.ULabel.filter(name__in=names).count() == 0
