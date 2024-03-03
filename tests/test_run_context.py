import lamindb as ln
from lamindb.core._run_context import get_uid_ext, run_context


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
    ln.core.run_context.run = None
    ln.core.run_context.transform = None


def test_track_with_reference():
    transform = ln.Transform(name="test")
    ln.track(transform, reference="my address", reference_type="url")
    assert ln.core.run_context.run.reference == "my address"
    assert ln.core.run_context.run.reference_type == "url"
    # unset to remove side effects
    ln.core.run_context.run = None
    ln.core.run_context.transform = None


def test_track_notebook_colab():
    notebook_path = "/fileId=1KskciVXleoTeS_OGoJasXZJreDU9La_l"
    ln.core.run_context._track_notebook(path=notebook_path)


def test_create_or_load_transform():
    title = "title"
    stem_uid = "NJvdsWWbJlZS"
    version = "0"
    uid = "NJvdsWWbJlZS6K79"
    assert uid == f"{stem_uid}{get_uid_ext(version)}"
    run_context._create_or_load_transform(
        stem_uid=stem_uid,
        version=version,
        name=title,
        transform_type="notebook",
    )
    assert run_context.transform.uid == uid
    assert run_context.transform.version == version
    run_context._create_or_load_transform(
        transform=run_context.transform,
        stem_uid=stem_uid,
        version=version,
        name=title,
    )
    assert run_context.transform.uid == uid
    assert run_context.transform.version == version
