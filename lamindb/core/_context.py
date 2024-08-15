from __future__ import annotations

import builtins
import hashlib
import os
from datetime import datetime, timezone
from pathlib import Path, PurePath
from typing import TYPE_CHECKING

from lamin_utils import logger
from lamindb_setup.core.hashing import hash_file
from lnschema_core import Run, Transform, ids
from lnschema_core.ids import base62_12
from lnschema_core.users import current_user_id

from ._settings import settings
from ._sync_git import get_transform_reference_from_git_repo
from ._track_environment import track_environment
from .exceptions import (
    MissingContext,
    NotebookNotSavedError,
    NoTitleError,
    UpdateContext,
)
from .subsettings._transform_settings import transform_settings
from .versioning import bump_version as bump_version_function
from .versioning import increment_base62

if TYPE_CHECKING:
    from lamindb_setup.core.types import UPathStr
    from lnschema_core.types import TransformType

is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)

msg_path_failed = (
    "failed to infer notebook path.\nfix: pass `path` to ln.context.track()"
)


def get_uid_ext(version: str) -> str:
    from lamin_utils._base62 import encodebytes

    # merely zero-padding the nbproject version such that the base62 encoding is
    # at least 4 characters long doesn't yields sufficiently diverse hashes and
    # leads to collisions; it'd be nice because the uid_ext would be ordered
    return encodebytes(hashlib.md5(version.encode()).digest())[:4]  # noqa: S324


def update_stem_uid_or_version(
    stem_uid: str | None,
    version: str | None,
    bump_version: bool | None = None,
) -> (bool, str, str):  # type:ignore
    updated = False
    if bump_version:
        response = "bump"
    elif bump_version is None:
        # ask for generating a new transform
        # it simply looks better here to not use the logger because we won't have an
        # emoji also for the subsequent input question
        if os.getenv("LAMIN_TESTING") is None:
            response = input(
                "To create a new transform, type 'new'. To bump the version, type 'bump'"
                " or a custom version: "
            )
        else:
            response = "new"
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
            new_version = bump_version_function(version, behavior="prompt")
        else:
            new_version = response
        updated = new_version != version
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
        name = get(f"http://{ip}:9000/api/sessions").json()[0]["name"]  # noqa: S113
    except Exception:
        logger.warning(
            "could not get notebook name from Google Colab, using: notebook.ipynb"
        )
        name = "notebook.ipynb"
    return name.rstrip(".ipynb")


def raise_missing_context(transform_type: str, key: str) -> None:
    transform = Transform.filter(key=key).latest_version().first()
    if transform is None:
        new_uid = f"{base62_12()}0000"
        message = f"To track this {transform_type}, set"
    else:
        uid = transform.uid
        suid, vuid = uid[: Transform._len_stem_uid], uid[Transform._len_stem_uid :]
        new_vuid = increment_base62(vuid)
        new_uid = f"{suid}{new_vuid}"
        suffix = ".ipynb" if transform_type == "notebook" else ".py"
        message = f"You already have a {transform_type} with filename '{key}{suffix}'. To create a new {transform_type}, rename your file. To bump the version, set"
    message += f'\n\nln.context.uid = "{new_uid}"'
    raise MissingContext(message)


def pretty_pypackages(dependencies: dict) -> str:
    deps_list = []
    for pkg, ver in dependencies.items():
        if ver != "":
            deps_list.append(pkg + f"=={ver}")
        else:
            deps_list.append(pkg)
    deps_list.sort()
    return " ".join(deps_list)


class Context:
    """Run context."""

    def __init__(self):
        self._uid: str | None = None
        self._name: str | None = None
        self._version: str | None = None
        self._transform: Transform | None = None
        self._run: Run | None = None
        self._path: Path | None = None
        """A local path to the script that's running."""
        self._logging_message: str = ""

    @property
    def transform(self) -> Transform | None:
        """Transform of context."""
        if self._transform is not None:
            return self._transform
        elif self._run is not None:
            return self._run.transform
        else:
            return None

    @property
    def uid(self) -> str | None:
        """`uid` to create transform."""
        return self._uid

    @uid.setter
    def uid(self, value: str | None):
        self._uid = value

    @property
    def name(self) -> str | None:
        """`name` to create transform."""
        return self._name

    @name.setter
    def name(self, value: str | None):
        self._name = value

    @property
    def version(self) -> str | None:
        """`version` to create transform."""
        return self._version

    @version.setter
    def version(self, value: str | None):
        self._version = value

    @property
    def run(self) -> Run | None:
        """Run of context."""
        return self._run

    def track(
        self,
        *,
        params: dict | None = None,
        transform: Transform | None = None,
        new_run: bool | None = None,
        path: str | None = None,
    ) -> None:
        """Track notebook or script run.

        Creates or loads a global :class:`~lamindb.Run` that enables data
        lineage tracking. You can find it in :class:`~lamindb.core.context`.

        Saves source code and compute environment.

        If :attr:`~lamindb.core.Settings.sync_git_repo` is set, will first check
        whether the script exists in the git repository and add a link.

        Args:
            params: A dictionary of parameters to track for the run.
            transform: Can be of type `"pipeline"` or `"notebook"`
                (:class:`~lamindb.core.types.TransformType`).
            new_run: If `False`, loads latest run of transform
                (default notebook), if `True`, creates new run (default pipeline).
            path: Filepath of notebook or script. Only needed if it can't be
                automatically detected.

        Examples:

            To track a notebook or script, call:

            >>> import lamindb as ln
            >>> ln.context.track()

            If you'd like to track an abstract pipeline run, pass a
            :class:`~lamindb.Transform` object of ``type`` ``"pipeline"``:

            >>> ln.Transform(name="Cell Ranger", version="2", type="pipeline").save()
            >>> transform = ln.Transform.filter(name="Cell Ranger", version="2").one()
            >>> ln.context.track(transform=transform)
        """
        self._path = None
        if transform is None:
            is_tracked = False
            transform_settings_are_set = (
                transform_settings.stem_uid is not None
                and transform_settings.version is not None
            )
            transform = None
            stem_uid = None
            if self.uid is not None:
                transform = Transform.filter(uid=self.uid).one_or_none()
                if (
                    transform is not None
                    and transform.version is not None
                    and self.version is not None
                    and self.version != transform.version
                ):
                    logger.warning(
                        "loaded transform version doesn't match passed version"
                    )
            elif transform_settings_are_set:
                stem_uid, self.version = (
                    transform_settings.stem_uid,
                    transform_settings.version,
                )
                transform = Transform.filter(
                    uid__startswith=stem_uid, version=self.version
                ).one_or_none()
            if is_run_from_ipython:
                key, name = self._track_notebook(path=path)
                transform_type = "notebook"
                transform_ref = None
                transform_ref_type = None
            else:
                (name, key, transform_ref, transform_ref_type) = self._track_script(
                    path=path
                )
                transform_type = "script"
            if self.uid is not None or transform_settings_are_set:
                # overwrite whatever is auto-detected in the notebook or script
                if self.name is not None:
                    name = self.name
                self._create_or_load_transform(
                    uid=self.uid,
                    stem_uid=stem_uid,
                    version=self.version,
                    name=name,
                    transform_ref=transform_ref,
                    transform_ref_type=transform_ref_type,
                    transform_type=transform_type,
                    key=key,
                    transform=transform,
                )
                # if no error is raised, the transform is tracked
                is_tracked = True
            if not is_tracked:
                raise_missing_context(transform_type, key)
        else:
            if transform.type in {"notebook", "script"}:
                raise ValueError(
                    "Use ln.context.track() without passing transform in a notebook or script"
                    " - metadata is automatically parsed"
                )
            transform_exists = None
            if transform.id is not None:
                # transform has an id but unclear whether already saved
                transform_exists = Transform.filter(id=transform.id).first()
            if transform_exists is None:
                transform.save()
                self._logging_message += f"created Transform('{transform.uid}')"
                transform_exists = transform
            else:
                self._logging_message += f"loaded Transform('{transform.uid}')"
            self._transform = transform_exists

        if new_run is None:  # for notebooks, default to loading latest runs
            new_run = False if self._transform.type == "notebook" else True  # type: ignore

        run = None
        if not new_run:  # try loading latest run by same user
            run = (
                Run.filter(transform=self._transform, created_by_id=current_user_id())
                .order_by("-created_at")
                .first()
            )
            if run is not None:  # loaded latest run
                run.started_at = datetime.now(timezone.utc)  # update run time
                self._logging_message += f" & loaded Run('{run.started_at}')"

        if run is None:  # create new run
            run = Run(
                transform=self._transform,
                params=params,
            )
            self._logging_message += f" & created Run('{run.started_at}')"
        # can only determine at ln.finish() if run was consecutive in
        # interactive session, otherwise, is consecutive
        run.is_consecutive = True if is_run_from_ipython else None
        # need to save in all cases
        run.save()
        if params is not None:
            run.params.add_values(params)
        self._run = run
        track_environment(run)
        logger.important(self._logging_message)
        self._logging_message = ""

    def _track_script(
        self,
        *,
        path: UPathStr | None,
    ) -> tuple[str, str, str, str]:
        if path is None:
            import inspect

            frame = inspect.stack()[2]
            module = inspect.getmodule(frame[0])
            self._path = Path(module.__file__)
        else:
            self._path = Path(path)
        name = self._path.name
        key = name
        reference = None
        reference_type = None
        if settings.sync_git_repo is not None:
            reference = get_transform_reference_from_git_repo(self._path)
            reference_type = "url"
        return name, key, reference, reference_type

    def _track_notebook(
        self,
        *,
        path: str | None,
    ):
        if path is None:
            path = get_notebook_path()
        key = Path(path).stem
        if isinstance(path, (Path, PurePath)):
            path_str = path.as_posix()  # type: ignore
        else:
            path_str = str(path)
        if path_str.endswith("Untitled.ipynb"):
            raise RuntimeError("Please rename your notebook before tracking it")
        if path_str.startswith("/fileId="):
            key = get_notebook_name_colab()
            name = key
        else:
            import nbproject

            try:
                nbproject_title = nbproject.meta.live.title
            except IndexError:
                raise NotebookNotSavedError(
                    "The notebook is not saved, please save the notebook and"
                    " rerun `ln.context.track()`"
                ) from None
            if nbproject_title is None:
                raise NoTitleError(
                    "Please add a title to your notebook in a markdown cell: # Title"
                ) from None
            name = nbproject_title
        # log imported python packages
        if not path_str.startswith("/fileId="):
            try:
                from nbproject.dev._pypackage import infer_pypackages

                nb = nbproject.dev.read_notebook(path_str)
                logger.important(
                    "notebook imports:"
                    f" {pretty_pypackages(infer_pypackages(nb, pin_versions=True))}"
                )
            except Exception:
                logger.debug("inferring imported packages failed")
                pass
        self._path = Path(path_str)
        return key, name

    def _create_or_load_transform(
        self,
        *,
        uid: str | None,
        stem_uid: str | None,
        version: str | None,
        name: str,
        transform_ref: str | None = None,
        transform_ref_type: str | None = None,
        key: str | None = None,
        transform_type: TransformType = None,
        transform: Transform | None = None,
    ):
        # make a new transform record
        if transform is None:
            if uid is None:
                uid = f"{stem_uid}{get_uid_ext(version)}"
            transform = Transform(
                uid=uid,
                version=version,
                name=name,
                key=key,
                reference=transform_ref,
                reference_type=transform_ref_type,
                type=transform_type,
            )
            transform.save()
            self._logging_message += f"created Transform('{transform.uid}')"
        else:
            uid = transform.uid
            # check whether the transform file has been renamed
            if transform.key != key:
                new_suid = ids.base62_12()
                transform_type = "Notebook" if is_run_from_ipython else "Script"
                raise UpdateContext(
                    f"{transform_type} filename changed, create new transform by setting:\n"
                    f'ln.context.uid = "{new_suid}0000"'
                )
            elif transform.name != name:
                transform.name = name
                transform.save()
                self._logging_message += f"updated Transform('{transform.uid}')"
            # check whether transform source code was already saved
            if transform._source_code_artifact_id is not None:
                response = None
                if is_run_from_ipython:
                    response = "y"  # auto-bump version
                else:
                    hash, _ = hash_file(self._path)  # ignore hash_type for now
                    if hash != transform._source_code_artifact.hash:
                        response = "y"  # auto-bump version
                    else:
                        self._logging_message += f"loaded Transform('{transform.uid}')"
                if response is not None:
                    change_type = (
                        "Re-running saved notebook"
                        if is_run_from_ipython
                        else "Source code changed"
                    )
                    suid, vuid = (
                        uid[: Transform._len_stem_uid],
                        uid[Transform._len_stem_uid :],
                    )
                    new_vuid = increment_base62(vuid)
                    raise UpdateContext(
                        f"{change_type}, bump version by setting:\n\n"
                        f'ln.context.uid = "{suid}{new_vuid}"'
                    )
            else:
                self._logging_message += f"loaded Transform('{transform.uid}')"
        self._transform = transform


context = Context()
