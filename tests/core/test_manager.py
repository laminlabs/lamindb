import lamindb as ln


def test_manager_list():
    label = ln.ULabel(name="manager label")
    label.save()
    label_names = [f"ULabel {i}" for i in range(3)]
    labels = [ln.ULabel(name=name) for name in label_names]
    ln.save(labels)
    label.parents.set(labels)
    assert len(label.parents.to_list()) == 3
    assert "ULabel 1" in label.parents.to_list("name")
    label.delete(permanent=True)
    for label in labels:
        label.delete(permanent=True)
