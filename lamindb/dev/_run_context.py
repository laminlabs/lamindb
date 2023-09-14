import builtins
import hashlib
import re
from datetime import datetime, timezone
from pathlib import Path, PurePath
from typing import Dict, List, Optional, Tuple, Union

from lamin_utils import logger
from lamindb_setup import settings
from lamindb_setup.dev import InstanceSettings
from lnschema_core import Run, Transform
from lnschema_core.types import TransformType

from lamindb.dev.versioning import get_ids_from_old_version, init_id

from .hashing import to_b64_str

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


def save_or_load(msg):
    if logger._verbosity > 2:
        logger.success(msg)
    else:
        logger.important(
            msg.replace("loaded: ", "").replace("saved: ", "").replace("updated: ", "")
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


# also see lamindb-setup._notebook for related code
def reinitialize_notebook(
    transform: Transform, metadata: Optional[Dict] = None, ask_for_new_id: bool = True
) -> Tuple[Transform, Dict]:
    from nbproject import dev as nb_dev
    from nbproject._header import _filepath

    # init new_id, new_version
    new_id, new_version = transform.id, None

    if ask_for_new_id:
        response = input("Do you want to generate a new id? (y/n) ")
    else:
        response = "y"
    if response == "y":
        new_id = init_id(n_full_id=14)
    else:
        response = input("Do you want to set a new version (e.g. '1.1')? (y/n) ")
        if response == "y":
            new_version = input("Please type the version: ")
            new_id = get_ids_from_old_version(
                is_new_version_of=transform, version=new_version, n_full_id=14
            )

    nb = None
    if metadata is None:
        nb = nb_dev.read_notebook(_filepath)
        metadata = nb.metadata["nbproject"]

    # nbproject metadata needs the shorter nbproject id
    metadata["id"] = new_id[:12]
    metadata["version"] = new_version

    # here we check that responses to both inputs (for new id and version) were not 'n'
    if transform.id != new_id or transform.version != new_version:
        transform = Transform(
            id=new_id, version=new_version, type=TransformType.notebook
        )
    return transform, metadata


# from https://stackoverflow.com/questions/61901628
def get_notebook_name_colab() -> str:
    from socket import gethostbyname, gethostname  # type: ignore

    from requests import get  # type: ignore

    ip = gethostbyname(gethostname())  # 172.28.0.12
    try:
        name = get(f"http://{ip}:9000/api/sessions").json()[0]["name"]
    except Exception:
        logger.warning(
            "could not get notebook name from Google Colab, using: Notebook.ipynb"
        )
        name = "Notebook.ipynb"
    return name.rstrip(".ipynb")


def get_transform_kwargs_from_nbproject(
    nbproject_id: str, nbproject_version: str, nbproject_title: str
) -> Tuple[Optional[Transform], str, str, str, Optional[Transform]]:
    id_ext = to_b64_str(hashlib.md5(nbproject_version.encode()).digest())[:2]
    id = nbproject_id + id_ext
    version = nbproject_version
    transform = Transform.filter(
        id__startswith=nbproject_id, version=version
    ).one_or_none()
    name = nbproject_title
    old_version_of = None
    if transform is None:
        old_version_of = Transform.filter(id__startswith=nbproject_id).first()
    return transform, id, version, name, old_version_of


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
        reference: Optional[str] = None,
        reference_type: Optional[str] = None,
        notebook_path: Optional[str] = None,
        pypackage: Optional[Union[str, List[str]]] = None,
        editor: Optional[str] = None,
    ) -> None:
        """Track global `Transform` & `Run` for a notebook or pipeline.

        Creates or loads a :class:`~lamindb.Run` record and sets a global
        :class:`~lamindb.dev.run_context`.

        In a Jupyter notebook, call without any argument (metadata is parsed).
        If the notebook has no associated metadata ("is not initialized"),
        attempts to write metadata to disk. If it fails to so interactively, it
        will ask you to leverage the CLI.

        Args:
            transform: Can be of type `"pipeline"` or `"notebook"`
                (:class:`~lamindb.dev.types.TransformType`).
            new_run: If `False`, loads latest run of transform
                (default notebook), if True, creates new run (default pipeline).
            reference: Reference to pass to :class:`~lamindb.Run` record.
            reference_type: Reference type to pass to :class:`~lamindb.Run`
                record (e.g. "url").
            notebook_path: Filepath of notebook. Only needed if inference fails.
            pypackage: One or more python packages for which to parse versions.
            editor: Editor environment.
                Pass `'lab'` for jupyter lab and `'notebook'` for jupyter notebook,
                this can help to identify the correct mechanism for interactivity
                when automatic inference fails.

        Examples:

            If you're in a Jupyter notebook and installed lamindb with `pip
            install[jupyter]`, you can simply call:

            >>> ln.track()
            >>> ln.dev.run_context.transform
            >>> ln.dev.run_context.run

            If you'd like to track a pipeline, pass a
            :class:`~lamindb.Transform` object of `type` `"pipeline"`:

            >>> ln.Transform(name="Cell Ranger", version="2", type="pipeline").save()
            >>> transform = ln.Transform.filter(name="Cell Ranger", version="2").one()
            >>> ln.track(transform)
        """
        cls.instance = settings.instance
        import lamindb as ln

        if transform is None:
            is_tracked_notebook = False

            if is_run_from_ipython:
                try:
                    cls._track_notebook(
                        pypackage=pypackage,
                        notebook_path=notebook_path,
                        editor=editor,
                    )
                    is_tracked_notebook = True
                except Exception as e:
                    if isinstance(e, ImportError):
                        logger.warning(
                            "it looks like you are running ln.track() from a "
                            "notebook!\nplease install nbproject: pip install nbproject"
                        )
                    elif isinstance(e, UpdateNbWithNonInteractiveEditorError):
                        raise e
                    elif isinstance(e, (NotebookNotSavedError, NoTitleError)):
                        raise e
                    else:
                        logger.warning(f"automatic tracking of notebook failed: {e}")
                        raise e
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
                save_or_load(f"saved: {transform}")
                transform_exists = transform
            else:
                save_or_load(f"loaded: {transform}")
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
                run.reference = reference
                run.reference_type = reference_type
                run.save()
                save_or_load(f"loaded: {run}")

        if run is None:  # create new run
            run = ln.Run(
                transform=cls.transform,
                reference=reference,
                reference_type=reference_type,
            )
            run.save()
            save_or_load(f"saved: {run}")
        cls.run = run

        # at this point, we have a transform can display its parents if there are any
        parents = cls.transform.parents.all() if cls.transform is not None else []
        if len(parents) > 0:
            if len(parents) == 1:
                save_or_load(f"  parent transform: {parents[0]}")
            else:
                parents_formatted = "\n   - ".join([f"{parent}" for parent in parents])
                save_or_load(f"  parent transforms:\n   - {parents_formatted}")

        # only for newly intialized notebooks
        if hasattr(cls, "_notebook_meta"):
            _write_notebook_meta(cls._notebook_meta)  # type: ignore
            del cls._notebook_meta  # type: ignore

    @classmethod
    def _track_notebook(
        cls,
        *,
        notebook_path: Optional[str] = None,
        pypackage: Optional[Union[str, List[str]]] = None,
        editor: Optional[str] = None,
    ):
        """Infer Jupyter notebook metadata and create `Transform` record.

        Args:
            pypackage: One or more python packages to track.
            notebook_path: Filepath of notebook. Only needed if automatic
                inference fails.
            editor: Editor environment. Only needed if automatic inference fails.
                Pass `'lab'` for jupyter lab and `'notebook'` for jupyter notebook,
                this can help to identify the correct mechanism for interactivity
                when automatic inference fails.
        """
        import nbproject
        from nbproject.dev._jupyter_communicate import (
            notebook_path as get_notebook_path,
        )

        cls.instance = settings.instance

        metadata = None
        needs_init = False
        is_interactive = False
        reference = None
        version = None
        colab_id = None
        nbproject_id = None
        if notebook_path is None:
            path_env = None
            try:
                path_env = get_notebook_path(return_env=True)
            except Exception:
                raise RuntimeError(msg_path_failed)
            if path_env is None:
                raise RuntimeError(msg_path_failed)
            notebook_path, _env = path_env
        else:
            notebook_path, _env = notebook_path, editor
        if isinstance(notebook_path, (Path, PurePath)):
            notebook_path_str = notebook_path.as_posix()  # type: ignore
        else:
            notebook_path_str = str(notebook_path)
        if notebook_path_str.endswith("Untitled.ipynb"):
            raise RuntimeError("Please rename your notebook before tracking it")
        if notebook_path_str.startswith("/fileId="):
            # This is Google Colab!
            # google colab fileID looks like this
            # /fileId=1KskciVXleoTeS_OGoJasXZJreDU9La_l
            # we'll take the first 12 characters
            colab_id = notebook_path_str.replace("/fileId=", "")
            reference = f"colab_id: {colab_id}"
            filestem = get_notebook_name_colab()
            _env = "colab"
        else:
            try:
                metadata, needs_init, nb = nbproject.header(
                    pypackage=pypackage,
                    filepath=notebook_path_str,
                    env=_env if editor is None else editor,
                    metadata_only=True,
                )
                # this contains filepath if the header was run successfully
                from nbproject._header import _env, _filepath  # type: ignore
            except Exception as e:
                nbproject_failed_msg = (
                    "Auto-retrieval of notebook name & title failed.\n\nFixes: Either"
                    f" init on the CLI `lamin track {notebook_path_str}` or pass"
                    " transform manually `ln.track(ln.Transform(name='My"
                    " notebook'))`\n\nPlease consider pasting error at:"
                    f" https://github.com/laminlabs/nbproject/issues/new\n\n{e}"
                )
                raise RuntimeError(nbproject_failed_msg)
            try:
                from nbproject.dev._metadata_display import DisplayMeta
                from nbproject.dev._pypackage import infer_pypackages

                dm = DisplayMeta(metadata)
                logger.important(
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
                msg = msg_manual_init.format(notebook_path=notebook_path_str)
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

            nbproject_id = metadata["id"]
            nbproject_version = metadata["version"]
            filestem = Path(_filepath).stem
            try:
                nbproject_title = nbproject.meta.live.title
            except IndexError:
                raise NotebookNotSavedError(
                    "The notebook is not saved, please save the notebook and"
                    " rerun `ln.track()`"
                )
            if nbproject_title is None:
                raise NoTitleError(
                    "Please add a title to your notebook in a markdown cell: # Title"
                )
        # colab parsing succesful
        if colab_id is not None:
            id = colab_id[:14]
            transform = Transform.filter(id=id).one_or_none()
            name = filestem
            short_name = None
            old_version_of = None
        # nbproject parsing successful
        elif nbproject_id is not None:
            (
                transform,
                id,
                version,
                name,
                old_version_of,
            ) = get_transform_kwargs_from_nbproject(
                nbproject_id, nbproject_version, nbproject_title
            )
            short_name = filestem
        # make a new transform record
        if transform is None:
            transform = Transform(
                id=id,
                version=version,
                name=name,
                short_name=short_name,
                reference=reference,
                is_new_version_of=old_version_of,
                type=TransformType.notebook,
            )
            transform.save()
            save_or_load(f"saved: {transform}")
        else:
            # check whether there was an update
            if transform.name != name or transform.short_name != short_name:
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
                transform.name = name
                transform.short_name = filestem
                transform.save()
                if response == "y":
                    save_or_load(f"saved: {transform}")
                else:
                    save_or_load(f"updated: {transform}")
            else:
                save_or_load(f"loaded: {transform}")
        cls.transform = transform
