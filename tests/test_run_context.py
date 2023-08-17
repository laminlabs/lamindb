import lamindb as ln


def test_track_with_multi_parents():
    parent1 = ln.Transform(name="parent 1")
    parent1.save()
    parent2 = ln.Transform(name="parent 2")
    parent2.save()
    child = ln.Transform(name="Child")
    child.save()
    child.parents.set([parent1, parent2])
    ln.track(child, reference="my address", reference_type="url")
    # unset to remove side effects
    ln.dev.run_context.run = None
    ln.dev.run_context.transform = None


def test_track_with_reference():
    transform = ln.Transform(name="test")
    ln.track(transform, reference="my address", reference_type="url")
    assert ln.dev.run_context.run.reference == "my address"
    assert ln.dev.run_context.run.reference_type == "url"
    # unset to remove side effects
    ln.dev.run_context.run = None
    ln.dev.run_context.transform = None
