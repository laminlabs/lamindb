import re
from pathlib import Path
from typing import List, Optional, Union

import lnschema_core
import nbproject
from lamin_logger import logger
from lndb import settings
from lndb.dev import InstanceSettings
from lnschema_core import Run, Transform
from nbproject._is_run_from_ipython import is_run_from_ipython

msg_init_complete = (
    "Init complete. Hit save & reload from disk, i.e, *discard* editor content. If you"
    " do not want to lose editor changes, hit save *before* running `track()`."
    " Consider using Jupyter Lab for a seamless interactive experience."
)


def _write_notebook_meta(metadata):
    from nbproject._header import _env, _filepath

    nbproject.dev._frontend_commands._save_notebook(_env)
    nb = nbproject.dev.read_notebook(_filepath)
    nb.metadata["nbproject"] = metadata

    # write proper execution count
    header_re = re.compile(r"^[^#]*track\(", flags=re.MULTILINE)
    ccount = 0
    for cell in nb.cells:
        if cell["cell_type"] != "code":
            continue
        elif cell["execution_count"] is not None:
            ccount = cell["execution_count"]
        if header_re.match("".join(cell["source"])) is not None:
            cell["execution_count"] = ccount + 1
            break

    nbproject.dev.write_notebook(nb, _filepath)
    nbproject.dev._frontend_commands._reload_notebook(_env)


class context:
    """Global run context."""

    instance: Optional[InstanceSettings] = None
    """Current instance."""
    transform: Optional[Transform] = None
    """Current transform."""
    run: Optional[Run] = None
    """Current run."""

    @classmethod
    def _track(
        cls,
        *,
        transform: Optional[Transform] = None,
        load_latest_run: bool = False,
        notebook_path: Optional[str] = None,
        pypackage: Union[str, List[str], None] = None,
        editor: Optional[str] = None,
    ) -> None:
        """Track `Transform` & `Run` records for a notebook or pipeline.

        Adds these records to the DB and exposes them as
        `ln.context.transform` and `ln.context.run`.

        Call without a `transform` record or without arguments
        when tracking a Jupyter notebook.

        If a Jupyter notebook has no associated metadata, attempts to write
        metadata to disk.

        Args:
            transform: Can be "pipeline" or "notebook".
            load_latest_run: If True, loads latest run of transform.
            pypackage: One or more python packages to track.
            notebook_path: Filepath of notebook. Only needed if inference fails.
            editor: Editor environment. Only needed if automatic inference fails.
                Pass `'lab'` for jupyter lab and `'notebook'` for jupyter notebook,
                this can help to identify the correct mechanism for interactivity
                when automatic inference fails.
        """
        cls.instance = settings.instance
        logger.info(f"Instance: {cls.instance.identifier}")
        logger.info(f"User: {settings.user.handle}")
        import lamindb as ln

        if is_run_from_ipython and transform is None:
            cls._track_notebook(
                pypackage=pypackage,
                filepath=notebook_path,
                editor=editor,
            )
        elif transform is None:
            raise ValueError("Pass `transform` to .track()!")
        else:
            transform_exists = ln.select(Transform, id=transform.id).one_or_none()
            if transform_exists is None:
                transform_exists = ln.add(transform)
                logger.info(f"Added transform: {transform}")
            else:
                logger.info(f"Loaded transform: {transform_exists}")
            cls.transform = transform_exists

        # this here uses cls.transform and writes cls.run
        # should probably change that design
        Run(load_latest=load_latest_run)

        # only for newly intialized notebooks
        if hasattr(cls, "_notebook_meta"):
            _write_notebook_meta(cls._notebook_meta)  # type: ignore
            del cls._notebook_meta  # type: ignore

    @classmethod
    def _track_notebook(
        cls,
        *,
        id: Optional[str] = None,
        v: Optional[str] = "0",
        name: Optional[str] = None,
        filepath: Optional[str] = None,
        pypackage: Optional[Union[str, List[str]]] = None,
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
        # original location of this code was _nb
        # legacy code here, see duplicated version in _run
        if id is None and name is None:
            nbproject_failed_msg = (
                "Auto-retrieval of notebook name & title failed.\nPlease paste error"
                " at: https://github.com/laminlabs/nbproject/issues/new \n\nFix: Run"
                " ln.track(transform=ln.Transform(name='My notebook',"
                " type='notebook'))"
            )
            try:
                metadata, needs_init = nbproject.header(
                    pypackage=pypackage,
                    filepath=filepath,
                    env=editor,
                    metadata_only=True,
                )
            except Exception:
                raise RuntimeError(nbproject_failed_msg)
            # this contains filepath if the header was run successfully
            from nbproject._header import _env, _filepath

            if needs_init:
                if _env in ("lab", "notebook"):
                    cls._notebook_meta = metadata  # type: ignore
                else:
                    nb = nbproject.dev.read_notebook(_filepath)
                    nb.metadata["nbproject"] = metadata
                    nbproject.dev.write_notebook(nb, _filepath)
                    raise SystemExit(msg_init_complete)

            id = metadata["id"]
            v = metadata["version"]
            name = Path(_filepath).stem
            title = nbproject.meta.live.title
        elif id is None or name is None:
            raise RuntimeError(
                "Cannot infer Transform metadata from notebook, pass `transform` to"
                " ln.track()."
            )
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
    def _track_pipeline(
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
                    f"Did not find any transform record with name '{name}'. Create a"
                    " new one? (y/n)"
                )
                if response == "y":
                    transform = Transform(name=name, type="pipeline")
        cls.transform = transform
