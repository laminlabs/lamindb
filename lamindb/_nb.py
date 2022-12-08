from typing import List, Union

import nbproject as nb
from lamin_logger import logger
from lnschema_core import Notebook, Run

_notebook: Notebook = None  # Notebook of this Python session
_run: Run = None  # run of this Python session


def header(
    *,
    run: Union[str, None] = None,
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

    notebook = ln.select(
        lns.Notebook, id=nb.meta.store.id, v=nb.meta.store.version
    ).one_or_none()
    if notebook is None:
        notebook = lns.Notebook(
            id=nb.meta.store.id, v=nb.meta.store.version, name=nb.meta.live.title
        )
        notebook = ln.add(notebook)
        logger.info(f"Added notebook: {notebook.id} v{notebook.v}")

    # at this point, we have a notebook object
    global _notebook
    _notebook = notebook

    # check user input
    # if isinstance(run, lns.Run):
    # This here might be something we may want in the future
    # but checking all the cases in which that run record has integrity
    # is quite a bit of code - not now!
    #     run_test = ln.select(lns.Run, id=run.id).one_or_none()
    #     if run_test is None:
    #         logger.info("Passed run does not exist, adding it")
    #         ln.add(run)
    if run is None:
        # retrieve the latest run
        run = (
            ln.select(lns.Run, notebook_id=notebook.id, notebook_v=notebook.v)
            .order_by(lns.Run.created_at.desc())
            .first()
        )
        if run is not None:
            logger.info(f"Loaded run: {run.id}")  # type: ignore
    elif run != "new":
        raise ValueError("Pass run='new' to header()!")

    # create a new run if doesn't exist yet or is requested by the user ("new")
    if run is None or run == "new":
        run = lns.Run(notebook_id=notebook.id, notebook_v=notebook.v)
        run = ln.add(run)
        logger.info(f"Created run: {run.id}")  # type: ignore

    # at this point, we have a run object
    global _run
    _run = run


def publish(version: str = None, i_confirm_i_saved: bool = False):
    """Publish the notebook.

    Args:
        version: Notebook version to publish.
        i_confirm_i_saved: Only relevant outside Jupyter Lab & Notebook as a
            safeguard against losing the editor buffer content.
    """
    from nbproject._publish import finalize_publish, run_checks_for_publish

    import lamindb as ln
    import lamindb.schema as lns

    result = run_checks_for_publish(
        calling_statement="publish(", i_confirm_i_saved=i_confirm_i_saved
    )
    if result != "checks-passed":
        return result
    finalize_publish(calling_statement="publish(", version=version)
    # update DB
    notebook = ln.select(Notebook, id=_notebook.id, v=_notebook.v).one()
    # update version
    notebook.name = nb.meta.live.title
    if version != _notebook.v:
        notebook_add = lns.Notebook(id=notebook.id, v=version, name=notebook.name)
    else:
        notebook_add = notebook
    ln.add(notebook_add)
    if version != _notebook.v:
        _run.notebook_v = version
        ln.add(_run)
        ln.delete(notebook)
