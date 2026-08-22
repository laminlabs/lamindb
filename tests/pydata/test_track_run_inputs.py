import lamindb as ln
import pytest


@pytest.fixture
def create_dataset():
    """Factory fixture that returns a function to create artifacts and collections."""
    created_sqlrecords = []

    def create(kind: str) -> ln.models.SQLRecord:
        if kind == "artifact":
            sqlrecord = ln.Artifact("README.md", key="README.md").save()
        elif kind == "collection":
            a1 = ln.Artifact("README.md", key="README.md").save()
            created_sqlrecords.append(a1)
            a2 = ln.Artifact("pyproject.toml", key="pyproject.toml").save()
            created_sqlrecords.append(a2)
            sqlrecord = ln.Collection([a1, a2], key="test-collection").save()
        created_sqlrecords.append(sqlrecord)
        return sqlrecord

    yield create

    for sqlrecord in created_sqlrecords[::-1]:
        sqlrecord.delete(permanent=True)


@pytest.mark.parametrize("registry_str", ["artifact", "collection"])
def test_track_run_input(create_dataset, registry_str):
    # First run - create the artifact/collection
    ln.track()
    # create an object
    sqlrecord = create_dataset(registry_str)
    # .cache() triggers input tracking
    sqlrecord.cache()
    # here tracking this object as in input of the current is skipped
    # because it was just created and we would get a cycle between the input and the output
    assert sqlrecord not in getattr(ln.context.run, f"input_{registry_str}s").all()
    # store the current global run, we will need it later
    first_run = ln.context.run
    first_sqlrecord = sqlrecord

    # Second run -- recreate the artifact/collection
    ln.track()
    # the new global run is not the same as the previous one
    assert ln.context.run != first_run
    # create a new artifact or collection, which will trigger a hash look up
    # and return the same sqlrecord as before
    sqlrecord = create_dataset(registry_str)
    assert sqlrecord == first_sqlrecord
    # because that run created the same sqlrecord, it is tracked as a recreating run
    assert ln.context.run in sqlrecord.recreating_runs.all()
    # we also track it in the private attribute to avoid database queries
    assert sqlrecord._recreating_run_id == ln.context.run.id
    # when we now trigger input tracking it's actually skipped
    # because we would create a cycle between the input and the fact that this artifact/collection
    # was recreated in this run
    # the skipping mechanism is cheap because it works through the cached _recreating_run_id attribute
    sqlrecord.cache()
    # assert that there is indeed no cycle
    assert sqlrecord not in getattr(ln.context.run, f"input_{registry_str}s").all()

    # Third run - retrieve the artifact/collection
    ln.track()
    assert ln.context.run != first_run
    # now we're querying the object
    if registry_str == "artifact":
        sqlrecord = ln.Artifact.get(key="README.md")
    else:
        sqlrecord = ln.Collection.get(key="test-collection")
    # trigger input tracking by calling .cache()
    sqlrecord.cache()
    # now it's tracked that this sqlrecord is an input of the current run
    assert sqlrecord in getattr(ln.context.run, f"input_{registry_str}s").all()
    # this run does not re-create this sqlrecord
    assert ln.context.run not in sqlrecord.recreating_runs.all()
    assert not hasattr(sqlrecord, "_recreating_run_id")
    # attempt to re-create the sqlrecord after it was retrieved in the same run
    sqlrecord = create_dataset(registry_str)
    assert sqlrecord == first_sqlrecord
    # because that run already registered this sqlrecord as an input
    # it is still not tracked as a recreating run because
    # we'd otherwise create a cycle
    assert ln.context.run not in sqlrecord.recreating_runs.all()
    assert not hasattr(sqlrecord, "_recreating_run_id")
