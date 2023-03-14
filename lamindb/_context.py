from pathlib import Path
from typing import List, Optional, Union

import lnschema_core
import nbproject
from lamin_logger import logger
from lndb import settings
from lndb.dev import InstanceSettings
from lnschema_core import Notebook, Pipeline, Run, dev
from nbproject._is_run_from_ipython import is_run_from_ipython


class context:
    """Global run context.

    Set through `ln.Run(global_context=True)`.

    Often, you'll want to call: `ln.Run(global_context=True, load_latest)`.
    """

    instance: Optional[InstanceSettings] = None
    """Current instance."""
    notebook: Optional[Notebook] = None
    """Current notebook."""
    pipeline: Optional[Pipeline] = None
    """Current pipeline."""
    run: Optional[Run] = None
    """Current run."""

    @classmethod
    def _track_notebook(
        cls,
        *,
        id: Optional[str] = None,
        v: Optional[str] = "0",
        name: Optional[str] = None,
        filepath: Optional[str] = None,
        pypackage: Union[str, List[str], None] = None,
        editor: Optional[str] = None,
    ):
        """Track notebook.

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
        cls.instance = settings.instance
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
            logger.info(f"Added notebook: {notebook}")
        else:
            logger.info(f"Loaded notebook: {notebook}")
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
    def _track_pipeline(
        cls,
        name: str,
        *,
        version: Optional[str] = None,
    ):
        """Track pipeline.

        Args:
            name: Pipeline name.
            version: Pipeline version. If `None`, load latest (sort by created_at).
        """
        cls.instance = settings.instance
        import lamindb as ln
        import lamindb.schema as lns

        if version is not None:
            pipeline = ln.select(lns.Pipeline, name=name, v=version).one()
        else:
            pipeline = (
                ln.select(lns.Pipeline, name=name)
                .order_by(lns.Pipeline.created_at.desc())
                .first()
            )
            if pipeline is None:
                response = input(
                    f"Did not find any pipeline record with name '{name}'. Create a new"
                    " one? (y/n)"
                )
                if response == "y":
                    pipeline = lns.Pipeline(name=name)
        cls.pipeline = pipeline

    @classmethod
    def _track_notebook_pipeline(
        cls, *, pipeline_name: Optional[str] = None, load_latest=True
    ):
        """Track notebook/pipeline and run.

        When called from within a Python script, pass `pipeline_name`.

        Args:
            pipeline_name: Pipeline name.
            load_latest: Load the latest run of the notebook or pipeline.
        """
        cls.instance = settings.instance
        logger.info(f"Instance: {cls.instance.identifier}")
        logger.info(f"User: {settings.user.handle}")
        if is_run_from_ipython and pipeline_name is None:
            if context.notebook is None:
                cls._track_notebook()
        else:
            if pipeline_name is None:
                raise ValueError(
                    "Pass a pipeline name: ln.context.track(pipeline_name='...')"
                )
            cls._track_pipeline(name=pipeline_name)
            logger.info(f"Pipeline: {cls.pipeline}")
