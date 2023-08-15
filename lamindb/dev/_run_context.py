import builtins
import os
import re
from datetime import datetime, timezone
from pathlib import Path, PurePath
from typing import Dict, List, Optional, Tuple, Union

import lnschema_core
from lamin_utils import logger
from lamindb_setup import settings
from lamindb_setup.dev import InstanceSettings
from lnschema_core import Run, Transform
from lnschema_core.types import TransformType

is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)

msg_path_failed = (
    "failed to infer notebook path.\nfix: either track manually via"
    " `ln.track(ln.Transform(name='My notebook'))` or pass"
    " `notebook_path` to ln.track()"
)

msg_manual_init = (
    "\n(1) save your notebook!"
    "\n(2) attach metadata to the notebook by running the CLI:\n"
    "lamin track {notebook_path}"
    "\n(3) reload or re-open your notebook"
)


class UpdateNbWithNonInteractiveEditorError(Exception):
    pass


class NotebookNotSavedError(Exception):
    pass


class NoTitleError(Exception):
    pass


def _seconds_modified(filepath):
    return datetime.now().timestamp() - Path(filepath).stat().st_mtime


def _write_notebook_meta(metadata):
    from nbproject import dev as nb_dev
    from nbproject._header import _env, _filepath

    nb_dev._frontend_commands._save_notebook(_env)
    nb = nb_dev.read_notebook(_filepath)
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

    nb_dev.write_notebook(nb, _filepath)
    nb_dev._frontend_commands._reload_notebook(_env)


def reinitialize_notebook(
    transform: Transform, metadata: Optional[Dict] = None
) -> Tuple[Transform, Dict]:
    from nbproject import dev as nb_dev
    from nbproject._header import _filepath

    new_id, new_version = transform.stem_id, None

    if "NBPRJ_TEST_NBPATH" not in os.environ:
        response = input("Do you want to generate a new id? (y/n)")
    else:
        response = "y"
    if response == "y":
        new_id = lnschema_core.ids.base62_12()
    else:
        response = input("Do you want to set a new version (e.g. '1.1')? (y/n)")
        if response == "y":
            new_version = input("Please type the version: ")

    nb = None
    if metadata is None:
        nb = nb_dev.read_notebook(_filepath)
        metadata = nb.metadata["nbproject"]

    metadata["id"] = new_id
    if new_version is None:
        new_version = "0"
    metadata["version"] = new_version

    # here we check that responses to both inputs (for new id and version) were not 'n'
    if transform.stem_id != new_id or transform.version != new_version:
        transform = Transform(
            stem_id=new_id, version=new_version, type=TransformType.notebook
        )
    return transform, metadata


# from https://stackoverflow.com/questions/61901628
def get_notebook_name_colab() -> str:
    from socket import gethostbyname, gethostname  # type: ignore

    from requests import get  # type: ignore

    ip = gethostbyname(gethostname())  # 172.28.0.12
    name = get(f"http://{ip}:9000/api/sessions").json()[0]["name"]
    return name.rstrip(".ipynb")


class run_context:
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
        """Track global `Transform` & `Run` for a notebook or pipeline.

        Access them via `ln.context.transform` and `ln.context.run`.

        Call without a `transform` record or without arguments
        when tracking a Jupyter notebook.

        If a Jupyter notebook has no associated metadata, attempts to write
        metadata to disk.

        Args:
            transform: `Optional[Transform] = None` - Can be of type
                "pipeline" or "notebook" (:class:`lamindb.types.TransformType`).
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

        Examples:

            If you're in a Jupyter notebook and installed lamindb with `pip
            install[jupyter]`, you can simply call:

            >>> ln.track()
            ✅ saved: Transform(id=1LCd8kco9lZUBg, name=Track data lineage / provenance, short_name=02-data-lineage, stem_id=1LCd8kco9lZU, version=0, type=notebook, updated_at=2023-07-10 18:37:19, created_by_id=DzTjkKse) # noqa
            ✅ saved: Run(id=pHgVICV9DxBaV6BAuKJl, run_at=2023-07-10 18:37:19, transform_id=1LCd8kco9lZUBg, created_by_id=DzTjkKse) # noqa
            >>> ln.context.transform
            Transform(id=1LCd8kco9lZUBg, name=Track data lineage / provenance, short_name=02-data-lineage, stem_id=1LCd8kco9lZU, version=0, type=notebook, updated_at=2023-07-10 18:37:19, created_by_id=DzTjkKse) # noqa
            >>> ln.context.run
            Run(id=pHgVICV9DxBaV6BAuKJl, run_at=2023-07-10 18:37:19, transform_id=1LCd8kco9lZUBg, created_by_id=DzTjkKse) # noqa

            If you'd like to track a pipeline we need to pass a
            :class:`~lamindb.Transform` object of `type` `"pipeline"`:

            >>> ln.Transform(name="Cell Ranger", version="7.2.0", type="pipeline").save()
            >>> transform = ln.Transform.filter(name="Cell Ranger", version="7.2.0").one()
            >>> ln.track(transform)
            💬 loaded: Transform(id=ceHkZMaiHFdoB6, name=Cell Ranger, stem_id=ceHkZMaiHFdo, version=7.2.0, type=pipeline, updated_at=2023-07-10 18:37:19, created_by_id=DzTjkKse) # noqa
            ✅ saved: Run(id=RcpWIKC8cF74Pn3RUJ1W, run_at=2023-07-10 18:37:19, transform_id=ceHkZMaiHFdoB6, created_by_id=DzTjkKse) # noqa
            >>> ln.context.transform
            Transform(id=ceHkZMaiHFdoB6, name=Cell Ranger, stem_id=ceHkZMaiHFdo, version=7.2.0, type=pipeline, updated_at=2023-07-10 18:37:19, created_by_id=DzTjkKse) # noqa
            >>> ln.context.run
            Run(id=RcpWIKC8cF74Pn3RUJ1W, run_at=2023-07-10 18:37:19, transform_id=ceHkZMaiHFdoB6, created_by_id=DzTjkKse) # noqa
        """
        cls.instance = settings.instance
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
                            "it looks like you are running ln.track() from a "
                            "notebook!\nplease install nbproject: pip install nbproject"
                        )
                    elif isinstance(e, UpdateNbWithNonInteractiveEditorError):
                        raise e
                    elif isinstance(e, (NotebookNotSavedError, NoTitleError)):
                        raise e
                    else:
                        logger.warning(f"automatic tracking of notebook failed: {e}")
                    is_tracked_notebook = False

            if not is_tracked_notebook:
                logger.warning(
                    "no automatic metadata detection, consider passing transform"
                )
                return None
        else:
            transform_exists = None
            if transform.id is not None:
                # transform has an id but unclear whether already saved
                transform_exists = Transform.filter(id=transform.id).first()
            if transform_exists is None:
                transform.save()
                logger.save(f"saved: {transform}")
                transform_exists = transform
            else:
                logger.success(f"loaded: {transform_exists}")
            cls.transform = transform_exists

        if new_run is None:  # for notebooks, default to loading latest runs
            new_run = False if cls.transform.type == TransformType.notebook.value else True  # type: ignore  # noqa

        run = None
        if not new_run:  # try loading latest run by same user
            run = (
                ln.Run.filter(
                    transform=cls.transform, created_by_id=ln.setup.settings.user.id
                )
                .order_by("-created_at")
                .first()
            )
            if run is not None:  # loaded latest run
                run.run_at = datetime.now(timezone.utc)  # update run time
                run.save()
                logger.success(f"loaded: {run}")

        if run is None:  # create new run
            run = ln.Run(transform=cls.transform)
            run.save()
            logger.save(f"saved: {run}")
        cls.run = run

        # at this point, we have a transform can display its parents if there are any
        parents = cls.transform.parents.all() if cls.transform is not None else []
        if len(parents) > 0:
            if len(parents) == 1:
                logger.info(f"  parent transform: {parents[0]}")
            else:
                parents_formatted = "\n   - ".join([f"{parent}" for parent in parents])
                logger.info(f"  parent transforms:\n   - {parents_formatted}")

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
        is_interactive = False
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
            notebook_path, _env = filepath, editor
        if isinstance(notebook_path, (Path, PurePath)):
            notebook_path = notebook_path.as_posix()
        if notebook_path.endswith("Untitled.ipynb"):
            raise RuntimeError("Please rename your notebook before tracking it")
        if notebook_path.startswith("/fileId="):
            # This is Google Colab!
            # google colab fileID looks like this
            # /fileId=1KskciVXleoTeS_OGoJasXZJreDU9La_l
            # we'll take the first 12 characters
            colab_id = notebook_path.replace("/fileId=", "")
            id = colab_id[:12]
            reference = f"colab_id: {colab_id}"
            filestem = get_notebook_name_colab()
            _env = "colab"
        else:
            try:
                metadata, needs_init, nb = nbproject.header(
                    pypackage=pypackage,
                    filepath=notebook_path if filepath is None else filepath,
                    env=_env if editor is None else editor,
                    metadata_only=True,
                )
                # this contains filepath if the header was run successfully
                from nbproject._header import _env, _filepath  # type: ignore
            except Exception as e:
                nbproject_failed_msg = (
                    "Auto-retrieval of notebook name & title failed.\n\nFixes: Either"
                    f" init on the CLI `lamin track {notebook_path}` or pass"
                    " transform manually `ln.track(ln.Transform(name='My"
                    " notebook'))`\n\nPlease consider pasting error at:"
                    f" https://github.com/laminlabs/nbproject/issues/new\n\n{e}"
                )
                raise RuntimeError(nbproject_failed_msg)
            try:
                from nbproject.dev._metadata_display import DisplayMeta
                from nbproject.dev._pypackage import infer_pypackages

                dm = DisplayMeta(metadata)
                logger.info(
                    "notebook imports:"
                    f" {' '.join(dm.pypackage(infer_pypackages(nb, pin_versions=True)))}"  # noqa
                )
            except Exception:
                logger.debug("inferring imported packages failed")
                pass

        if needs_init:
            if _env in ("lab", "notebook"):
                cls._notebook_meta = metadata  # type: ignore
            else:
                msg = msg_manual_init.format(notebook_path=notebook_path)
                raise UpdateNbWithNonInteractiveEditorError(msg)

        if _env in ("lab", "notebook"):
            # save the notebook in case that title was updated
            # but notebook not saved
            nbproject.dev._frontend_commands._save_notebook(_env)
            # check here if interactivity really works

        if metadata is not None:
            # this only executed if nbproject is used
            # check here if interactivity really works
            # timestamp should be recent due to save
            is_interactive = _seconds_modified(_filepath) < 1.5  # should be ~1 sec
            if not is_interactive and needs_init:
                msg = msg_manual_init.format(notebook_path=_filepath)
                raise UpdateNbWithNonInteractiveEditorError(msg)

            id = metadata["id"]
            version = metadata["version"]
            filestem = Path(_filepath).stem
            try:
                title = nbproject.meta.live.title
            except IndexError:
                raise NotebookNotSavedError(
                    "The notebook is not saved, please save the notebook and"
                    " rerun `ln.track()`"
                )
            if title is None:
                raise NoTitleError(
                    "Please add a title to your notebook in a markdown cell: # Title"
                )
        else:
            version = "0"
            title = filestem

        transform = Transform.filter(stem_id=id, version=version).one_or_none()
        if transform is None:
            transform = Transform(
                stem_id=id,
                version=version,
                name=title,
                short_name=filestem,
                reference=reference,
                type=TransformType.notebook,
            )
            transform.save()
            logger.save(f"saved: {transform}")
        else:
            logger.success(f"loaded: {transform}")
            if transform.name != title or transform.short_name != filestem:
                response = input(
                    "Updated notebook name and/or title: Do you want to assign a"
                    " new id or version? (y/n)"
                )
                if response == "y":
                    if _env in ("lab", "notebook") and is_interactive:
                        transform, metadata = reinitialize_notebook(transform, metadata)
                        # only write metadata back to notebook if it actually changed!
                        # if filename or title changed, this does not merit a write!
                        # it's dangerous to write unnecessarily
                        cls._notebook_meta = metadata  # type: ignore
                    else:
                        msg = msg_manual_init.format(notebook_path=notebook_path)
                        raise UpdateNbWithNonInteractiveEditorError(msg)
                transform.name = title
                transform.short_name = filestem
                transform.save()
                if response == "y":
                    logger.save(f"saved: {transform}")
                else:
                    logger.success(f"updated: {transform}")

        cls.transform = transform
