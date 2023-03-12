from pathlib import Path
from typing import List, Optional, Union

import lnschema_core
import nbproject
from lamin_logger import logger
from lndb.dev import InstanceSettings
from lnschema_core import Notebook, Pipeline, Run, dev


class context:
    """Global run context.

    This is helpful if you don't want to ingest data from concurrent pipeline runs.
    """

    instance: InstanceSettings = None  # current instance
    notebook: Notebook = None  # current notebook
    pipeline: Pipeline = None  # current pipeline
    run: Run = None  # run of this Python session

    @classmethod
    def track_notebook(
        cls,
        *,
        id: Optional[str] = None,
        v: Optional[str] = "0",
        name: Optional[str] = None,
        filepath: Optional[str] = None,
        pypackage: Union[str, List[str], None] = None,
        editor: Optional[str] = None,
    ):
        """Track a Notebook global record.

        Args:
            id: Pass a notebook id manually.
            v: Pass a notebook version manually.
            name: Pass a notebook name manually.
            pypackage: One or more python packages to track.
            filepath: Filepath of notebook. Only needed if automatic inference fails.
            editor: Editor environment. Only needed if automatic inference fails.
                Pass `'lab'` for jupyter lab and `'notebook'` for jupyter notebook,
                this can help to identify the correct mechanism for interactivity
                when automatic inference fails.
        """
        # original location of this code was _nb
        # legacy code here, see duplicated version in _run
        if id is None and name is None:
            nbproject_failed_msg = (
                "Auto-retrieval of notebook name & title failed.\nPlease paste error"
                " at: https://github.com/laminlabs/nbproject/issues/new \n\nFix: Run"
                f" ln.nb.header(id={dev.id.notebook()}, name='my-notebook-name')"
            )
            try:
                nbproject.header(
                    pypackage=pypackage, filepath=filepath, env=editor, display=False
                )
            except Exception:
                raise RuntimeError(nbproject_failed_msg)
            # this contains filepath if the header was run successfully
            from nbproject._header import _filepath

            id = nbproject.meta.store.id
            v = nbproject.meta.store.version
            name = Path(_filepath).stem
            title = nbproject.meta.live.title
        elif id is None or name is None:
            # Both id and name need to be passed if passing it manually
            raise RuntimeError("Fix: Pass both id & name to ln.nb.header().")
        else:
            title = None

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
                response = input(
                    "Updated notebook name and/or title: Do you want to assign a new id"
                    " or version? (y/n)"
                )
                if response == "y":
                    print("Notebook metadata will be re-initialized.")
                    new_id, new_v = None, None
                    response = input("Do you want to generate a new id? (y/n)")
                    if response == "y":
                        new_id = lnschema_core.dev.id.notebook()
                    response = input(
                        "Do you want to set a new version (e.g. '1.1')? Type 'n' for"
                        " 'no'. (version/n)"
                    )
                    if new_v != "n":
                        if new_v == "y":
                            response = input("Please type the version: ")
                        new_v = response
                    if new_id is not None or new_v is not None:
                        nbproject.meta.store.id = new_id
                        nbproject.meta.store.version = new_v
                        nbproject.meta.store.write()
                        # at this point, depending on the editor, the process
                        # might crash that is OK as upon re-running, the
                        # notebook will have new metadata and will be registered
                        # in the db in case the python process does not exit, we
                        # need a new Notebook record
                        notebook = lns.Notebook(id=id, v=v)

                notebook.name = name
                notebook.title = title
                ln.add(notebook)

        # at this point, we have a notebook object
        cls.notebook = notebook

    @classmethod
    def track_run(
        cls,
        *,
        load_latest: bool = False,
        run: Optional[Run] = None,
    ):
        """Track a Notebook global record.

        Args:
            run: If `None`, either create new run or load latest.
            load_latest: Loads the latest run of the notebook or pipeline.
        """
        import lamindb as ln
        import lamindb.schema as lns

        # check user input
        # if isinstance(run, lns.Run):
        # This here might be something we may want in the future
        # but checking all the cases in which that run record has integrity
        # is quite a bit of code - not now!
        #     run_test = ln.select(lns.Run, id=run.id).one_or_none()
        #     if run_test is None:
        #         logger.info("Passed run does not exist, adding it")
        #         ln.add(run)
        if run is None and load_latest:
            if cls.notebook is not None:
                select_stmt = ln.select(
                    lns.Run, notebook_id=cls.notebook.id, notebook_v=cls.notebook.v
                )
            elif cls.pipeline is not None:
                select_stmt = ln.select(
                    lns.Run, pipeline_id=cls.pipeline.id, pipeline_v=cls.pipeline.v
                )
            else:
                raise RuntimeError(
                    "Please call context.track_notebook or context.track_pipeline."
                )
            run = select_stmt.order_by(lns.Run.created_at.desc()).first()
            if run is not None:
                logger.info(f"Loaded run: {run.id}")  # type: ignore
            else:
                logger.info("Did not find any latest run. Creating new run.")

        # create a new run if doesn't exist yet or is requested by the user
        if run is None:
            if cls.notebook is not None:
                run = lns.Run(notebook_id=cls.notebook.id, notebook_v=cls.notebook.v)
            elif cls.pipeline is not None:
                run = lns.Run(pipeline_id=cls.pipeline.id, pipeline_v=cls.pipeline.v)
            else:
                raise RuntimeError(
                    "Please call context.track_notebook or context.track_pipeline."
                )
            run = ln.add(run)  # type: ignore
            logger.info(f"Added run: {run.id}")  # type: ignore

        # at this point, we have a run object
        cls.run = run
