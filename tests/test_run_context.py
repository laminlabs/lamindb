import lamindb as ln


def test_track_with_multi_parents():
    parent1 = ln.Transform(name="parent 1")
    parent1.save()
    parent2 = ln.Transform(name="parent 2")
    parent2.save()
    child = ln.Transform(name="Child")
    child.save()
    child.parents.set([parent1, parent2])
    ln.track(child)
    # unset to remove side effects
    ln.dev.run_context.run = None
    ln.dev.run_context.transform = None
