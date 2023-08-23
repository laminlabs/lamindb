import lamindb as ln


def test_delete():
    labels = [ln.Label(name=name) for name in ["label1", "label2", "label3"]]
    ln.save(labels)
    ln.delete(labels[0])
    ln.delete([labels[1], labels[2]])
    assert ln.Label.filter(name__in=["label1", "label2", "label3"]).count() == 0
