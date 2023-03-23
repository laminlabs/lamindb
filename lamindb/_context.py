from pathlib import Path
from typing import List, Optional, Union

import lnschema_core
import nbproject
from lamin_logger import logger
from lndb import settings
from lndb.dev import InstanceSettings
from lnschema_core import Run, Transform, dev


class context:
    """Global run context.

    If you set `ln.context.transform = transform`, the `transform` record is
    recognized whenever you create a run with `ln.Run(global_context=True)`.

    There are two convenience methods for creating `transform` records for
    Jupyter notebooks and pipelines, each.

    You can load the latest run if it doesn't exist via:
    `ln.Run(global_context=True, load_latest=True)`.
    """

    instance: Optional[InstanceSettings] = None
    """Current instance."""
    transform: Optional[Transform] = None
    """Current transform."""
    run: Optional[Run] = None
    """Current run."""

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
        """Infer Jupyter notebook metadata and create `Transform` record.

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
        logger.info(f"Instance: {cls.instance.identifier}")
        logger.info(f"User: {settings.user.handle}")
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

        transform = ln.select(
            Transform,
            id=id,
            v=v,
        ).one_or_none()
        if transform is None:
            transform = Transform(
                id=id,
                v=v,
                name=name,
                title=title,
                type="notebook",
            )
            transform = ln.add(transform)
            logger.info(f"Added notebook: {transform}")
        else:
            logger.info(f"Loaded notebook: {transform}")
            if transform.name != name or transform.title != title:
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
                        transform = Transform(id=id, v=v, type="notebook")

                transform.name = name
                transform.title = title
                ln.add(transform)

        # at this point, we have a transform object
        cls.transform = transform

    @classmethod
    def track_pipeline(
        cls,
        name: str,
        *,
        version: Optional[str] = None,
    ):
        """Load or create pipeline record within `Transform`.

        Args:
            name: Name as used in `Transform.name`.
            version: Pipeline version. If `None`, load latest (sort by `created_at`).
        """
        cls.instance = settings.instance
        import lamindb as ln

        if version is not None:
            transform = ln.select(Transform, name=name, v=version).one()
        else:
            transform = (
                ln.select(Transform, name=name)
                .order_by(Transform.created_at.desc())
                .first()
            )
            if transform is None:
                response = input(
                    f"Did not find any pipeline record with name '{name}'. Create a new"
                    " one? (y/n)"
                )
                if response == "y":
                    transform = Transform(name=name, type="pipeline")
        cls.transform = transform
