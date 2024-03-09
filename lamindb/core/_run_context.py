import builtins
import hashlib
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path, PurePath
from typing import Any, Dict, List, Optional, Tuple, Union

from lamin_utils import logger
from lamindb_setup import settings
from lamindb_setup.core import InstanceSettings
from lnschema_core import Run, Transform, ids
from lnschema_core.types import TransformType
from lnschema_core.users import current_user_id

from lamindb.core._transform_settings import transform_settings

is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)

msg_path_failed = (
    "failed to infer notebook path.\nfix: either track manually via"
    " `ln.track(ln.Transform(name='My notebook'))` or pass"
    " `path` to ln.track()"
)


class NotebookNotSavedError(Exception):
    pass


class NoTitleError(Exception):
    pass


def get_uid_ext(version: str) -> str:
    from lamin_utils._base62 import encodebytes

    # merely zero-padding the nbproject version such that the base62 encoding is
    # at least 4 characters long doesn't yields sufficiently diverse hashes and
    # leads to collisions; it'd be nice because the uid_ext would be ordered
    return encodebytes(hashlib.md5(version.encode()).digest())[:4]


def get_stem_uid_and_version_from_file(file_path: str) -> Tuple[str, str]:
    # line-by-line matching might be faster, but let's go with this for now
    with open(file_path) as file:
        content = file.read()

    if file_path.endswith(".py"):
        stem_uid_pattern = re.compile(
            r'\.transform\.stem_uid\s*=\s*["\']([^"\']+)["\']'
        )
        version_pattern = re.compile(r'\.transform\.version\s*=\s*["\']([^"\']+)["\']')
    elif file_path.endswith(".ipynb"):
        stem_uid_pattern = re.compile(
            r'\.transform\.stem_uid\s*=\s*\\["\']([^"\']+)\\["\']'
        )
        version_pattern = re.compile(
            r'\.transform\.version\s*=\s*\\["\']([^"\']+)\\["\']'
        )
    else:
        raise ValueError("Only .py and .ipynb files are supported.")

    # Search for matches in the entire file content
    stem_uid_match = stem_uid_pattern.search(content)
    version_match = version_pattern.search(content)

    # Extract values if matches are found
    stem_uid = stem_uid_match.group(1) if stem_uid_match else None
    version = version_match.group(1) if version_match else None

    if stem_uid is None or version is None:
        raise SystemExit(
            f"ln.transform.stem_uid and ln.transform.version aren't set in {file_path}\n"
            "Call ln.track() and copy/paste the output into the notebook"
        )
    return stem_uid, version


def update_stem_uid_or_version(
    stem_uid: str,
    version: str,
    bump_version: bool = False,
) -> (bool, str, str):  # type:ignore
    get_uid_ext(version)
    updated = False
    if bump_version:
        response = "bump"
    else:
        # ask for generating a new stem uid
        # it simply looks better here to not use the logger because we won't have an
        # emoji also for the subsequent input question
        if os.getenv("LAMIN_TESTING") is None:
            response = input(
                "To create a new stem uid, type 'new'. To bump the version, type 'bump'"
                " or a custom version: "
            )
        else:
            response = "new"
        if response == "new":
            new_stem_uid = ids.base62_12()
            updated = True
        else:
            bump_version = True
    new_version = version
    if bump_version:
        new_stem_uid = stem_uid
        if response == "bump":
            try:
                new_version = str(int(version) + 1)
            except ValueError:
                new_version = input(
                    f"The current version is '{version}' - please type the new"
                    " version: "
                )
        else:
            new_version = response
        updated = new_version != version
    if updated:
        new_metadata = (
            f'ln.transform.stem_uid = "{new_stem_uid}"\nln.transform.version ='
            f' "{new_version}"\n'
        )
        raise SystemExit(
            f"Please update your transform settings as follows:\n{new_metadata}"
        )
    return updated, new_stem_uid, new_version


def get_notebook_path():
    from nbproject.dev._jupyter_communicate import (
        notebook_path as get_notebook_path,
    )

    path = None
    try:
        path = get_notebook_path()
    except Exception:
        raise RuntimeError(msg_path_failed) from None
    if path is None:
        raise RuntimeError(msg_path_failed) from None
    return path


# from https://stackoverflow.com/questions/61901628
def get_notebook_name_colab() -> str:
    from socket import gethostbyname, gethostname  # type: ignore

    from requests import get  # type: ignore

    ip = gethostbyname(gethostname())  # 172.28.0.12
    try:
        name = get(f"http://{ip}:9000/api/sessions").json()[0]["name"]
    except Exception:
        logger.warning(
            "could not get notebook name from Google Colab, using: notebook.ipynb"
        )
        name = "notebook.ipynb"
    return name.rstrip(".ipynb")


MESSAGE = """To track this {transform_type}, set the following two global variables:

ln.transform.stem_uid = "{stem_uid}"
ln.transform.version = "{version}"
"""


def raise_transform_settings_error() -> None:
    from lnschema_core.ids import base62_12

    transform_type = "notebook" if is_run_from_ipython else "script"
    stem_uid = base62_12()
    version = "1"

    # backward compat: use the nbproject_id
    if is_run_from_ipython:
        from nbproject.dev import (
            MetaContainer,
            MetaStore,
            read_notebook,
        )
        from nbproject.dev._jupyter_communicate import (
            notebook_path as get_notebook_path,
        )

        filepath = get_notebook_path()
        nb = read_notebook(filepath)  # type: ignore
        nb_meta = nb.metadata
        if "nbproject" in nb_meta:
            meta_container = MetaContainer(**nb_meta["nbproject"])
            meta_store = MetaStore(meta_container, filepath)
            stem_uid, version = meta_store.id, meta_store.version
    raise SystemExit(
        MESSAGE.format(
            transform_type=transform_type, stem_uid=stem_uid, version=version
        )
    )


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
        path: Optional[str] = None,
    ) -> None:
        """Track global `Transform` & `Run` for a notebook or pipeline.

        Creates or loads a :class:`~lamindb.Run` record and sets a global
        :class:`~lamindb.core.run_context`.

        In a Jupyter notebook, call without any argument (metadata is parsed).
        If the notebook has no associated metadata ("is not initialized"),
        attempts to write metadata to disk. If it fails to so interactively, it
        will ask you to leverage the CLI.

        Args:
            transform: Can be of type `"pipeline"` or `"notebook"`
                (:class:`~lamindb.core.types.TransformType`).
            new_run: If `False`, loads latest run of transform
                (default notebook), if True, creates new run (default pipeline).
            reference: Reference to pass to :class:`~lamindb.Run` record.
            reference_type: Reference type to pass to :class:`~lamindb.Run`
                record (e.g. "url").
            path: Filepath of notebook or script. Only needed if it can't be
                automatically detected.

        Examples:

            To track a notebook or script, call:

            >>> import lamindb as ln
            >>> ln.track()
            # if global transform settings are not yet defined, this will ask you to set them
            # if they are defined, this will log the transform and its run

            If you'd like to track a pipeline run, pass a
            :class:`~lamindb.Transform` object of `type` `"pipeline"`:

            >>> ln.Transform(name="Cell Ranger", version="2", type="pipeline").save()
            >>> transform = ln.Transform.filter(name="Cell Ranger", version="2").one()
            >>> ln.track(transform)
        """
        cls.instance = settings.instance
        if transform is None:
            is_tracked = False
            transform_settings_are_set = (
                transform_settings.stem_uid is not None
                and transform_settings.version is not None
            )
            if transform_settings_are_set:
                stem_uid, version = (
                    transform_settings.stem_uid,
                    transform_settings.version,
                )
                transform = Transform.filter(
                    uid__startswith=stem_uid, version=version
                ).one_or_none()
                if is_run_from_ipython:
                    short_name, name, _ = cls._track_notebook(path=path)
                else:
                    import inspect

                    frame = inspect.stack()[1]
                    module = inspect.getmodule(frame[0])
                    name = Path(module.__file__).name  # type: ignore
                    short_name = name
                transform_type = (
                    TransformType.notebook
                    if is_run_from_ipython
                    else TransformType.script
                )
                cls._create_or_load_transform(
                    stem_uid=stem_uid,
                    version=version,
                    name=name,
                    reference=reference,
                    transform_type=transform_type,
                    short_name=short_name,
                    transform=transform,
                )
                # if no error is raised, the transform is tracked
                is_tracked = True
            if not is_tracked:
                raise_transform_settings_error()
        else:
            transform_exists = None
            if transform.id is not None:
                # transform has an id but unclear whether already saved
                transform_exists = Transform.filter(id=transform.id).first()
            if transform_exists is None:
                transform.save()
                logger.important(f"saved: {transform}")
                transform_exists = transform
            else:
                logger.important(f"loaded: {transform}")
            cls.transform = transform_exists

        if new_run is None:  # for notebooks, default to loading latest runs
            new_run = (
                False if cls.transform.type == TransformType.notebook.value else True
            )  # type: ignore

        run = None
        from lamindb._run import Run

        if not new_run:  # try loading latest run by same user
            run = (
                Run.filter(transform=cls.transform, created_by_id=current_user_id())
                .order_by("-created_at")
                .first()
            )
            if run is not None:  # loaded latest run
                run.run_at = datetime.now(timezone.utc)  # update run time
                run.reference = reference
                run.reference_type = reference_type
                run.save()
                logger.important(f"loaded: {run}")

        if run is None:  # create new run
            run = Run(
                transform=cls.transform,
                reference=reference,
                reference_type=reference_type,
            )
            run.save()
            logger.important(f"saved: {run}")
        cls.run = run

        from ._track_environment import track_environment

        track_environment(run)

        # at this point, we have a transform can display its parents if there are any
        parents = cls.transform.parents.all() if cls.transform is not None else []
        if len(parents) > 0:
            if len(parents) == 1:
                logger.info(f"  parent transform: {parents[0]}")
            else:
                parents_formatted = "\n   - ".join([f"{parent}" for parent in parents])
                logger.info(f"  parent transforms:\n   - {parents_formatted}")

    @classmethod
    def _track_notebook(
        cls,
        *,
        path: Optional[str],
    ):
        if path is None:
            path = get_notebook_path()
        short_name = Path(path).stem
        if isinstance(path, (Path, PurePath)):
            path_str = path.as_posix()  # type: ignore
        else:
            path_str = str(path)
        if path_str.endswith("Untitled.ipynb"):
            raise RuntimeError("Please rename your notebook before tracking it")
        if path_str.startswith("/fileId="):
            short_name = get_notebook_name_colab()
            name = short_name
        else:
            import nbproject

            try:
                nbproject_title = nbproject.meta.live.title
            except IndexError:
                raise NotebookNotSavedError(
                    "The notebook is not saved, please save the notebook and"
                    " rerun `ln.track()`"
                ) from None
            if nbproject_title is None:
                raise NoTitleError(
                    "Please add a title to your notebook in a markdown cell: # Title"
                ) from None
            name = nbproject_title
        # log imported python packages
        if not path_str.startswith("/fileId="):
            try:
                from nbproject.dev._metadata_display import DisplayMeta
                from nbproject.dev._pypackage import infer_pypackages

                metadata, _, nb = nbproject.header(
                    filepath=path_str,
                    metadata_only=True,
                )
                dm = DisplayMeta(metadata)
                logger.important(
                    "notebook imports:"
                    f" {' '.join(dm.pypackage(infer_pypackages(nb, pin_versions=True)))}"
                )
            except Exception:
                logger.debug("inferring imported packages failed")
                pass
        return short_name, name, path_str

    @classmethod
    def _create_or_load_transform(
        cls,
        *,
        stem_uid: str,
        version: Optional[str],
        name: str,
        reference: Optional[str] = None,
        short_name: Optional[str] = None,
        transform_type: TransformType = None,
        transform: Optional[Transform] = None,
    ):
        # make a new transform record
        if transform is None:
            uid = f"{stem_uid}{get_uid_ext(version)}"
            transform = Transform(
                uid=uid,
                version=version,
                name=name,
                short_name=short_name,
                reference=reference,
                type=transform_type,
            )
            transform.save()
            logger.important(f"saved: {transform}")
        else:
            # check whether there was an update to the transform, like
            # renaming the file or updating the title
            if transform.name != name or transform.short_name != short_name:
                if os.getenv("LAMIN_TESTING") is None:
                    response = input(
                        "Updated transform filename and/or title: Do you want to assign a"
                        " new stem_uid or version? (y/n)"
                    )
                else:
                    response = "y"
                if response == "y":
                    # will raise SystemExit
                    update_stem_uid_or_version(stem_uid, version)
                else:
                    transform.name = name
                    transform.short_name = short_name
                    transform.save()
                    logger.important(f"updated: {transform}")
            # check whether the transform artifacts were already saved
            if (
                transform.source_code_id is not None
                or transform.latest_report_id is not None
            ):
                if os.getenv("LAMIN_TESTING") is None:
                    response = input(
                        "You already saved a source file for this transform."
                        " Do you want to bump the version? (y/n)"
                    )
                else:
                    response = "y"
                if response == "y":
                    update_stem_uid_or_version(stem_uid, version, bump_version=True)
                else:
                    # we want a new stem_uid in this case, hence raise the error
                    raise_transform_settings_error()
            else:
                logger.important(f"loaded: {transform}")
        cls.transform = transform
