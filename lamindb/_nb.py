from typing import List, Union

import nbproject as nb
from lamin_logger import logger
from lnschema_core import Run
from nbproject import publish  # noqa

_run: Run = None  # run of this Python session


def header(
    *,
    run: Union[str, Run] = None,
    pypackage: Union[str, List[str], None] = None,
    filepath: Union[str, None] = None,
    env: Union[str, None] = None,
):
    """Display metadata and start tracking dependencies.

    If the notebook has no nbproject metadata, initializes & writes metadata to disk.

    Args:
        run: If `None`, loads the latest run of the notebook (or creates on if
            there isn't any).

            If `"new"`, instantiates a new `Run` record.

            If a `Run` record, loads that record.

        pypackage: One or more python packages to track.
        filepath: Filepath of notebook. Only needed if automatic inference fails.
        env: Editor environment. Only needed if automatic inference fails.
            Pass `'lab'` for jupyter lab and `'notebook'` for jupyter notebook,
            this can help to identify the correct mechanism for interactivity
            when automatic inference fails.
    """
    nb.header(pypackage=pypackage, filepath=filepath, env=env)

    import lamindb as ln
    import lamindb.schema as lns

    jupynb = ln.select(
        lns.Jupynb, id=nb.meta.store.id, v=nb.meta.store.version
    ).one_or_none()
    if jupynb is None:
        jupynb = lns.Jupynb(
            id=nb.meta.store.id, v=nb.meta.store.version, name=nb.meta.live.title
        )
        jupynb = ln.add(jupynb)
        logger.info(f"Added {jupynb}")

    # check user input
    if isinstance(run, lns.Run):
        run_test = ln.select(lns.Run, id=run.id).one_or_none()
        if run_test is None:
            logger.info("Passed run does not exist, adding it")
            ln.add(run)
    elif run is None:
        run = ln.select(lns.Run, jupynb_id=jupynb.id, jupynb_v=jupynb.v).one_or_none()
        if run is not None:
            logger.info(f"Loaded run: {run.id}")  # type: ignore
    elif run != "new":
        raise ValueError("Pass a lns.Run object to header() or 'new'!")

    # create a new run if doesn't exist yet or is requested by the user ("new")
    if run is None or run == "new":
        run = lns.Run(jupynb_id=jupynb.id, jupynb_v=jupynb.v)
        run = ln.add(run)
        logger.info(f"Added run: {run.id}")

    # at this point, we have a run object
    global _run
    _run = run
