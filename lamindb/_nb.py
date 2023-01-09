from pathlib import Path
from typing import List, Optional, Union

import nbproject as _nb
from lamin_logger import logger
from lndb_setup import settings
from lnschema_core import Notebook, Run, dev


class nb:
    """Manage Jupyter notebooks.

    For more background, see `nbproject <https://lamin.ai/docs/nbproject>`__.
    """

    notebook: Notebook = None  # Notebook of this Python session
    run: Run = None  # run of this Python session

    @classmethod
    def header(
        cls,
        *,
        run: Optional[str] = None,
        pypackage: Union[str, List[str], None] = None,
        filepath: Optional[str] = None,
        env: Optional[str] = None,
        id: Optional[str] = None,
        v: Optional[str] = "0",
        name: Optional[str] = None,
    ):
        """Track the notebook & display metadata.

        Call without arguments in most settings.

        If the notebook has no nbproject metadata, initializes & writes metadata
        to disk.

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
            id: Pass a notebook id manually.
            v: Pass a notebook version manually.
            name: Pass a notebook name manually.
        """
        if id is None and name is None:
            try:
                _nb.header(pypackage=pypackage, filepath=filepath, env=env)
            except Exception:
                raise RuntimeError(
                    "Failed to run nbproject.header(). "
                    f"Pass id={dev.id.notebook()} & name."
                )
            from nbproject.dev._jupyter_communicate import notebook_path

            id = _nb.meta.store.id
            v = _nb.meta.store.version
            name = Path(notebook_path()).stem
            title = _nb.meta.live.title
        elif id is None or name is None:
            # Both id and name need to be passed if passing it manually
            raise RuntimeError("Pass both id & name.")
        else:
            title = None

        logger.info(f"Instance: {settings.instance.owner}/{settings.instance.name}")

        import lamindb as ln
        import lamindb.schema as lns

        notebook = ln.select(
            lns.Notebook,
            id=id,
            v=v,
        ).one_or_none()
        if notebook is None:
            notebook = lns.Notebook(
                id=id,
                v=v,
                name=name,
                title=title,
            )
            notebook = ln.add(notebook)
            logger.info(f"Added notebook: {notebook.id} v{notebook.v}")
        else:
            logger.info(f"Loaded notebook: {notebook.id} v{notebook.v}")
            if notebook.name != name or notebook.title != title:
                notebook.name = name
                notebook.title = title
                ln.add(notebook)
                logger.info("Updated notebook name or title.")

        # at this point, we have a notebook object
        cls.notebook = notebook

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
            run = ln.add(run)  # type: ignore
            logger.info(f"Added run: {run.id}")  # type: ignore

        # at this point, we have a run object
        cls.run = run

    @classmethod
    def publish(cls, version: str = None, i_confirm_i_saved: bool = False):
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
        notebook = ln.select(Notebook, id=cls.notebook.id, v=cls.notebook.v).one()
        # update version
        notebook.title = _nb.meta.live.title
        if version != cls.notebook.v:
            notebook_add = lns.Notebook(
                id=notebook.id, v=version, name=notebook.name, title=notebook.title
            )
        else:
            notebook_add = notebook
        ln.add(notebook_add)
        if version != cls.notebook.v:
            cls.run.notebook_v = version
            ln.add(cls.run)
            ln.delete(notebook)
