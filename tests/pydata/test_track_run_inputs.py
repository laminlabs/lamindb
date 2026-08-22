import lamindb as ln
import pytest
from django.db.models import ProtectedError


@pytest.fixture
def create_dataset():
    """Factory fixture that returns a function to create artifacts and collections."""
    created_datasets = []

    def create(kind: str) -> ln.models.SQLRecord:
        if kind == "artifact":
            dataset = ln.Artifact("README.md", key="README.md").save()
        elif kind == "collection":
            a1 = ln.Artifact("README.md", key="README.md").save()
            created_datasets.append(a1)
            a2 = ln.Artifact("pyproject.toml", key="pyproject.toml").save()
            created_datasets.append(a2)
            dataset = ln.Collection([a1, a2], key="test-collection").save()
        created_datasets.append(dataset)
        return dataset

    yield create

    for dataset in created_datasets[::-1]:
        run_id = dataset.run_id
        recreating_runs = dataset.recreating_runs.all()
        for recreating_run in recreating_runs:
            recreating_run.delete(permanent=True)
        dataset.delete(permanent=True)
        if ln.Run.filter(id=run_id).exists():
            try:
                ln.Run.get(id=run_id).delete(permanent=True)
            except ProtectedError:
                pass
    ln.context._run = None


@pytest.mark.parametrize("registry_str", ["artifact", "collection"])
def test_track_run_input(create_dataset, registry_str):
    # First run - create the dataset
    ln.track()
    # create an object
    dataset = create_dataset(registry_str)
    # .cache() triggers input tracking
    dataset.cache()
    # here tracking this object as in input of the current is skipped
    # because it was just created and we would get a cycle between the input and the output
    assert dataset not in getattr(ln.context.run, f"input_{registry_str}s").all()
    # store the current global run, we will need it later
    first_run = ln.context.run
    first_dataset = dataset

    # Second run -- recreate the dataset
    ln.track()
    # the new global run is not the same as the previous one
    assert ln.context.run != first_run
    # create a new artifact or collection, which will trigger a hash look up
    # and return the same dataset as before
    dataset = create_dataset(registry_str)
    assert dataset == first_dataset
    # because that run created the same dataset, it is tracked as a recreating run
    assert ln.context.run in dataset.recreating_runs.all()
    # we also track it in the private attribute to avoid database queries
    assert dataset._recreating_run_id == ln.context.run.id
    # when we now trigger input tracking it's actually skipped
    # because we would create a cycle between the input and the fact that this artifact/collection
    # was recreated in this run
    # the skipping mechanism is cheap because it works through the cached _recreating_run_id attribute
    dataset.cache()
    # assert that there is indeed no cycle
    assert dataset not in getattr(ln.context.run, f"input_{registry_str}s").all()

    # Third run - retrieve the dataset
    ln.track()
    assert ln.context.run != first_run
    # now we're querying the object
    if registry_str == "artifact":
        dataset = ln.Artifact.get(key="README.md")
    else:
        dataset = ln.Collection.get(key="test-collection")
    # trigger input tracking by calling .cache()
    dataset.cache()
    # now it's tracked that this dataset is an input of the current run
    assert dataset in getattr(ln.context.run, f"input_{registry_str}s").all()
    # this run does not re-create this dataset
    assert ln.context.run not in dataset.recreating_runs.all()
    assert not hasattr(dataset, "_recreating_run_id")
    # attempt to re-create the dataset after it was retrieved in the same run
    dataset = create_dataset(registry_str)
    assert dataset == first_dataset
    # because that run already registered this dataset as an input
    # it is still not tracked as a recreating run because
    # we'd otherwise create a cycle
    assert ln.context.run not in dataset.recreating_runs.all()
    assert not hasattr(dataset, "_recreating_run_id")
