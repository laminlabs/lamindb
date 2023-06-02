import builtins
import os
import re
from pathlib import Path, PurePath
from typing import Dict, List, Optional, Tuple, Union

import lnschema_core
from lamin_logger import logger
from lamindb_setup import settings
from lamindb_setup.dev import InstanceSettings
from lnschema_core import Run, Transform

from .dev.db._add import _private_not_empty

is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)

msg_init_complete = (
    "⚠️ Destructive operation! ⚠️\n\nAre you sure you saved the notebook before running"
    " `ln.track()`?\n - If not, hit save, and *overwrite* the notebook file.\n - If"
    " yes, hit save, and *discard* editor content.\n\nConsider using Jupyter Lab or"
    " Notebook for a seamless interactive notebook tracking experience."
)

msg_path_failed = (
    "Failed to infer notebook path.\nFix: Either track manually via"
    " `ln.track(ln.Transform(name='My notebook'))` or pass"
    " `notebook_path` to ln.track()."
)


def _write_notebook_meta(metadata):
    import nbproject
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


def reinitialize_notebook(
    id: str, name: str, metadata: Optional[Dict] = None
) -> Tuple[Transform, Dict]:
    import nbproject
    from nbproject._header import _env, _filepath

    new_id, new_version = id, None
    if "NBPRJ_TEST_NBPATH" not in os.environ:
        response = input("Do you want to generate a new id? (y/n)")
    else:
        response = "y"
    if response == "y":
        new_id = lnschema_core.dev.id.transform()
    else:
        response = input(
            "Do you want to set a new version (e.g. '1.1')? Type 'n' for"
            " 'no'. (version/n)"
        )
        if response != "n":
            if new_version == "y":
                response = input("Please type the version: ")
            new_version = response

    nb = None
    if metadata is None:
        nb = nbproject.dev.read_notebook(_filepath)
        metadata = nb.metadata["nbproject"]

    metadata["id"] = new_id
    if new_version is None:
        new_version = "0"
    metadata["version"] = new_version

    # in "lab" & "notebook", we push the metadata write to the end of track execution
    # by returning metadata below
    if _env not in ("lab", "notebook", "test"):
        if nb is None:
            nb = nbproject.dev.read_notebook(_filepath)
        nb.metadata["nbproject"] = metadata
        nbproject.dev.write_notebook(nb, _filepath)
        raise SystemExit(msg_init_complete)

    transform = Transform(id=new_id, version=new_version, name=name, type="notebook")
    return transform, metadata


# from https://stackoverflow.com/questions/61901628
def get_notebook_name_colab() -> str:
    from socket import gethostbyname, gethostname  # type: ignore

    from requests import get  # type: ignore

    ip = gethostbyname(gethostname())  # 172.28.0.12
    name = get(f"http://{ip}:9000/api/sessions").json()[0]["name"]
    return name.rstrip(".ipynb")


class context:
    """Global run context."""

    instance: Optional[InstanceSettings] = None
    """Current instance."""
    transform: Optional[Transform] = None
    """Current transform."""
    run: Optional[Run] = None
    """Current run."""

    # exposed to user as ln.track()
    @classmethod
    def _track(
        cls,
        transform: Optional[Transform] = None,
        *,
        new_run: Optional[bool] = None,
        notebook_path: Optional[str] = None,
        pypackage: Optional[Union[str, List[str]]] = None,
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
            transform: `Optional[Transform] = None` - Can be of type
                "pipeline" or "notebook".
            new_run: `Optional[bool] = None` - If False, loads latest run of transform
                (default notebook), if True, creates new run (default pipeline).
            notebook_path: `Optional[str] = None` - Filepath of notebook.
                Only needed if inference fails.
            pypackage: `Optional[Union[str, List[str]]] = None` - One or more
                python packages to track.
            editor: `Optional[str] = None` - Editor environment.
                Pass `'lab'` for jupyter lab and `'notebook'` for jupyter notebook,
                this can help to identify the correct mechanism for interactivity
                when automatic inference fails.
        """
        cls.instance = settings.instance
        logger.info(f"Instance: {cls.instance.identifier}")
        logger.info(f"User: {settings.user.handle}")
        import lamindb as ln

        if transform is None:
            is_tracked_notebook = False

            if is_run_from_ipython:
                try:
                    cls._track_notebook(
                        pypackage=pypackage,
                        filepath=notebook_path,
                        editor=editor,
                    )
                    is_tracked_notebook = True
                except Exception as e:
                    if isinstance(e, ImportError):
                        logger.info(
                            "It looks like you are running ln.track() from a jupyter"
                            " notebook. Consider installing nbproject for automatic"
                            " tracking."
                        )
                    elif isinstance(e, RuntimeError):
                        raise e
                    is_tracked_notebook = False

            if not is_tracked_notebook:
                new_transform = Transform(name="Default pipeline", type="pipeline")
                ln.add(new_transform)
                logger.success(f"Added: {new_transform}")
                cls.transform = new_transform
        else:
            if transform.id is not None:  # id based look-up
                if transform.version is None:
                    transform_exists = (
                        ln.select(Transform, id=transform.id)
                        .order_by(Transform.created_at.desc())
                        .first()
                    )
                else:
                    transform_exists = ln.select(
                        Transform, id=transform.id, version=transform.version
                    ).first()
            else:  # name based lookup
                if transform.version is None:
                    transform_exists = (
                        ln.select(Transform, name=transform.name)
                        .order_by(Transform.created_at.desc())
                        .first()
                    )
                else:
                    transform_exists = ln.select(
                        Transform, name=transform.name, version=transform.version
                    ).first()
            if transform_exists is None:
                transform_exists = ln.add(transform)
                logger.success(f"Added: {transform}")
            else:
                logger.info(f"Loaded: {transform_exists}")
            cls.transform = transform_exists

        # for notebooks, default to loading latest runs
        if new_run is None:
            if cls.transform.type.value == "notebook":  # type: ignore
                new_run = False
            else:
                new_run = True

        # this here uses cls.transform and writes cls.run
        # should probably change that design
        Run(load_latest=not new_run)
        # so, this is a hack:
        if (
            _private_not_empty(cls.run, "_ln_identity_key")  # type: ignore
            and cls.run._ln_identity_key is not None  # type: ignore
        ):
            cls.run._sa_instance_state.key = cls.run._ln_identity_key  # type: ignore

        # only for newly intialized notebooks
        if hasattr(cls, "_notebook_meta"):
            _write_notebook_meta(cls._notebook_meta)  # type: ignore
            del cls._notebook_meta  # type: ignore

    @classmethod
    def _track_notebook(
        cls,
        *,
        filepath: Optional[str] = None,
        pypackage: Optional[Union[str, List[str]]] = None,
        editor: Optional[str] = None,
    ):
        """Infer Jupyter notebook metadata and create `Transform` record.

        Args:
            pypackage: One or more python packages to track.
            filepath: Filepath of notebook. Only needed if automatic inference fails.
            editor: Editor environment. Only needed if automatic inference fails.
                Pass `'lab'` for jupyter lab and `'notebook'` for jupyter notebook,
                this can help to identify the correct mechanism for interactivity
                when automatic inference fails.
        """
        import nbproject
        from nbproject.dev._jupyter_communicate import notebook_path

        cls.instance = settings.instance

        metadata = None
        needs_init = False
        reference = None
        if filepath is None:
            path_env = None
            try:
                path_env = notebook_path(return_env=True)
            except Exception:
                raise RuntimeError(msg_path_failed)
            if path_env is None:
                raise RuntimeError(msg_path_failed)
            notebook_path, _env = path_env
        else:
            notebook_path = filepath
        if isinstance(notebook_path, (Path, PurePath)):
            notebook_path = notebook_path.as_posix()
        if notebook_path.endswith("Untitled.ipynb"):
            raise RuntimeError("Please rename your notebook before tracking it")
        if notebook_path.startswith("/filedId="):
            # This is Google Colab!
            # google colab fileID looks like this
            # /fileId=1KskciVXleoTeS_OGoJasXZJreDU9La_l
            # we'll take the first 12 characters
            colab_id = notebook_path.replace("/filedId=", "")
            id = colab_id[:12]
            reference = f"colab_id: {colab_id}"
            name = get_notebook_name_colab()
            _env = "colab"
        else:
            try:
                metadata, needs_init = nbproject.header(
                    pypackage=pypackage,
                    filepath=notebook_path if filepath is None else filepath,
                    env=_env if editor is None else editor,
                    metadata_only=True,
                )
                # this contains filepath if the header was run successfully
                from nbproject._header import _env, _filepath  # type: ignore
            except Exception:
                nbproject_failed_msg = (
                    "Auto-retrieval of notebook name & title failed.\nPlease paste"
                    " error at: https://github.com/laminlabs/nbproject/issues/new"
                    " \n\nFix: Run `ln.track(ln.Transform(name='My notebook'))`"
                )
                raise RuntimeError(nbproject_failed_msg)

        import lamindb as ln

        if needs_init:
            if _env in ("lab", "notebook"):
                cls._notebook_meta = metadata  # type: ignore
            else:
                # nb = nbproject.dev.read_notebook(_filepath)
                # nb.metadata["nbproject"] = metadata
                # nbproject.dev.write_notebook(nb, _filepath)
                # raise SystemExit(msg_init_complete)
                # the following is safer
                raise RuntimeError(
                    "Please attach an ID to the notebook by running the CLI: lamin"
                    " track my-notebook.ipynb"
                )

        if _env in ("lab", "notebook"):
            # save the notebook in case that title was updated
            # but notebook not saved
            nbproject.dev._frontend_commands._save_notebook(_env)

        if metadata is not None:
            id = metadata["id"]
            version = metadata["version"]
            name = Path(_filepath).stem
            title = nbproject.meta.live.title
        else:
            version = "0"
            title = None

        transform = ln.select(Transform, id=id, version=version).one_or_none()
        if transform is None:
            transform = Transform(
                id=id,
                version=version,
                name=name,
                title=title,
                reference=reference,
                type="notebook",
            )
            transform = ln.add(transform)
            logger.success(f"Added: {transform}")
        else:
            logger.info(f"Loaded: {transform}")
            if transform.name != name or transform.title != title:
                if _env in ("lab", "notebook"):
                    response = input(
                        "Updated notebook name and/or title: Do you want to assign a"
                        " new id or version? (y/n)"
                    )
                    if response == "y":
                        transform, metadata = reinitialize_notebook(
                            transform.id, name, metadata
                        )
                    cls._notebook_meta = metadata  # type: ignore
                    transform.name = name
                    transform.title = title
                    ln.add(transform)
                    if response == "y":
                        logger.success(f"Added: {transform}")
                    else:
                        logger.success(f"Updated: {transform}")
                else:
                    logger.warning(
                        "Updated notebook name and/or title. If you want to assign a"
                        " new id or version, run: lamin track my-notebook.ipynb"
                    )

        cls.transform = transform
