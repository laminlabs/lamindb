# %% [markdown]
# # Loading relationships: `Session`
import pytest
from sqlalchemy.orm.exc import DetachedInstanceError

import lamindb as ln


def test_session():
    if ln._USE_DJANGO:
        return None

    # %% [markdown]
    # Let's create related sample data records and add them to the database:

    # %%
    transform = ln.Transform(name="Transform A")
    run = ln.Run(name="Solve Problem X", transform=transform)

    # %%
    transform

    # %%
    run

    # %%
    run.transform

    # %%
    ln.add(run)

    # %% [markdown]
    # Both records got just added to the database.

    # %% [markdown]
    # In the background, a `Session` object was created, which connected to the database, inserted the records, and closed the connection.

    # %% [markdown]
    # ## Query results without session

    # %%
    run_queried = ln.select(ln.Run, name="Solve Problem X").first()

    # %% [markdown]
    # Also here, in the background, a session was created and closed. This is good enough if we need to use simple properties of the returned record, for instance, the pipeline id:

    # %%
    run_queried

    # %% [markdown]
    # However, if we'd like to access the entire related record, we'll get a `DetachedInstanceError` error (it would tell us that the "lazy load operation of attribute 'inputs' cannot proceed"):

    # %%
    with pytest.raises(DetachedInstanceError):
        run_queried.inputs

    # %% [markdown]
    # The queried run would need to have an open connection to the DB in order for it to automatically load the related record. Under the hood, it needs to perform an automated query for this.
    #
    # But when `ln.select(...).first()` completed its execution, the database connection was closed.

    # %% [markdown]
    # ```{note}
    #
    # We can pre-configure to always load relationships in certain cases: `ln.Run.transform` is such a case, see [here](https://github.com/laminlabs/lnschema-core/blob/86a1984442ac514ada0aa3f246f2cab1022c795c/lnschema_core/_core.py#L161)!
    #
    # Hence, you'll be able to access `run_queried.transform`.
    #
    # ```
    #

    # %% [markdown]
    # ## The Session object

    # %% [markdown]
    # In order to lazily load related data records, we need to use a `Session` object!

    # %%
    ss = ln.Session()

    # %% [markdown]
    # The `Session` object comes with `add`, `delete` and `select`, just as the global namespace. They are equivalent to the global version, with the only difference being that all data records manipulated will be bound to an open session.

    # %%
    run_session = ss.select(ln.Run, name="Solve Problem X").first()

    # %%
    run_session

    # %% [markdown]
    # It's clear we don't need it for the simple attributes. But we need it for lazily loaded relationships:

    # %%
    run_session.inputs

    # %% [markdown]
    # Let us close the session.

    # %%
    ss.close()

    # %% [markdown]
    # Given we already loaded the pipeline record, it's still available in memory.

    # %%
    run_session.inputs

    # %% [markdown]
    # But, we can't access the `outputs` relationship, as the session is now closed.

    # %%
    with pytest.raises(DetachedInstanceError):
        run_session.outputs

    # %% [markdown]
    # ## The Session object in a context manager

    # %% [markdown]
    # We can also call `Session` in a context manager:

    # %%
    with ln.Session() as ss:
        run_session2 = ss.select(ln.Run, name="Solve Problem X").first()
        print(run_session2.outputs)

    # %% [markdown]
    # Because we loaded the ouputs, they're still in memory and available:

    # %%
    run_session2.outputs

    # %% [markdown]
    # Accessing another relationship, however, will error:

    # %%
    with pytest.raises(DetachedInstanceError):
        run_session2.inputs
