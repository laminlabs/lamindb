from __future__ import annotations

import builtins
import hashlib
from datetime import datetime, timezone
from pathlib import Path, PurePath
from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
from lamin_utils import logger
from lamindb_setup.core.hashing import hash_file
from lnschema_core import Run, Transform, ids
from lnschema_core.ids import base62_12
from lnschema_core.models import format_field_value

from ._settings import settings
from ._sync_git import get_transform_reference_from_git_repo
from ._track_environment import track_environment
from .exceptions import (
    InconsistentKey,
    MissingContextUID,
    NotebookNotSaved,
    NoTitleError,
    TrackNotCalled,
    UpdateContext,
)
from .subsettings._transform_settings import transform_settings
from .versioning import bump_version as bump_version_function
from .versioning import increment_base62, message_update_key_in_version_family

if TYPE_CHECKING:
    from lamindb_setup.core.types import UPathStr
    from lnschema_core.types import TransformType

is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)

msg_path_failed = "failed to infer notebook path.\nfix: pass `path` to `ln.track()`"


def get_uid_ext(version: str) -> str:
    from lamin_utils._base62 import encodebytes

    # merely zero-padding the nbproject version such that the base62 encoding is
    # at least 4 characters long doesn't yields sufficiently diverse hashes and
    # leads to collisions; it'd be nice because the uid_ext would be ordered
    return encodebytes(hashlib.md5(version.encode()).digest())[:4]  # noqa: S324


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


def raise_missing_context(transform_type: str, key: str) -> bool:
    transform = Transform.filter(key=key).latest_version().first()
    if transform is None:
        new_uid = f"{base62_12()}0000"
        message = f'to track this {transform_type}, copy & paste `ln.track("{new_uid}")` and re-run'
    else:
        uid = transform.uid
        new_uid = f"{uid[:-4]}{increment_base62(uid[-4:])}"
        message = f"you already have a transform with key '{key}' ('{transform.uid}')\n  - to make a revision, call `ln.track('{new_uid}')`\n  - to create a new transform, rename your file and run: `ln.track()`"
    if transform_type == "notebook":
        print(f"→ {message}")
        response = input("→ Ready to re-run? (y/n)")
        if response == "y":
            logger.important(
                "note: restart your notebook if you want consecutive cell execution"
            )
            return True
        raise MissingContextUID("Please follow the instructions.")
    else:
        raise MissingContextUID(f"✗ {message}")
    return False


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
    """Run context.

    Enables convenient data lineage tracking by managing a transform & run
    upon :meth:`~lamindb.core.Context.track` & :meth:`~lamindb.core.Context.finish`.

    Examples:

        Is typically used via the global :class:`~lamindb.context` object via `ln.track()` and `ln.finish()`:

        >>> import lamindb as ln
        >>> ln.track()
        >>> # do things
        >>> ln.finish()

    """

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
        """Managed transform of context."""
        return self._transform

    @property
    def uid(self) -> str | None:
        """`uid` argument for `context.transform`."""
        return self._uid

    @uid.setter
    def uid(self, value: str | None):
        self._uid = value

    @property
    def name(self) -> str | None:
        """`name argument for `context.transform`."""
        return self._name

    @name.setter
    def name(self, value: str | None):
        self._name = value

    @property
    def version(self) -> str | None:
        """`version` argument for `context.transform`."""
        return self._version

    @version.setter
    def version(self, value: str | None):
        self._version = value

    @property
    def run(self) -> Run | None:
        """Managed run of context."""
        return self._run

    def track(
        self,
        uid: str | None = None,
        *,
        params: dict | None = None,
        new_run: bool | None = None,
        path: str | None = None,
        transform: Transform | None = None,
    ) -> None:
        """Initiate a run with tracked data lineage.

        - sets :attr:`~lamindb.core.Context.transform` &
          :attr:`~lamindb.core.Context.run` by creating or loading `Transform` &
          `Run` records
        - saves compute environment as a `requirements.txt` file: `run.environment`

        If :attr:`~lamindb.core.Settings.sync_git_repo` is set, checks whether a
        script-like transform exists in a git repository and links it.

        Args:
            uid: A `uid` to create or load a transform.
            params: A dictionary of parameters to track for the run.
            new_run: If `False`, loads latest run of transform
                (default notebook), if `True`, creates new run (default pipeline).
            path: Filepath of notebook or script. Only needed if it can't be
                automatically detected.
            transform: Useful to track an abstract pipeline.

        Examples:

            To track the run of a notebook or script, call:

            >>> import lamindb as ln
            >>> ln.track()

        """
        if uid is not None:
            self.uid = uid
        self._path = None
        if transform is None:
            is_tracked = False
            transform_settings_are_set = (
                transform_settings.stem_uid is not None
                and transform_settings.version is not None
            )
            transform = None
            stem_uid = None
            if uid is not None or self.uid is not None:
                transform = Transform.filter(uid=self.uid).one_or_none()
                if self.version is not None:
                    # test inconsistent version passed
                    if (
                        transform is not None
                        and transform.version is not None
                        and self.version != transform.version
                    ):
                        raise SystemExit(
                            f"Please pass consistent version: ln.context.version = '{transform.version}'"
                        )
                    # test whether version was already used for another member of the family
                    suid, vuid = (
                        self.uid[: Transform._len_stem_uid],
                        self.uid[Transform._len_stem_uid :],
                    )
                    transform = Transform.filter(
                        uid__startswith=suid, version=self.version
                    ).one_or_none()
                    if (
                        transform is not None
                        and vuid != transform.uid[Transform._len_stem_uid :]
                    ):
                        better_version = bump_version_function(self.version)
                        raise SystemExit(
                            f"Version '{self.version}' is already taken by Transform(uid='{transform.uid}'); please set another version, e.g., ln.context.version = '{better_version}'"
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
                early_return = raise_missing_context(transform_type, key)
                if early_return:
                    return None
        else:
            if transform.type in {"notebook", "script"}:
                raise ValueError(
                    "Use `ln.track()` without passing transform in a notebook or script"
                    " - metadata is automatically parsed"
                )
            transform_exists = None
            if transform.id is not None:
                # transform has an id but unclear whether already saved
                transform_exists = Transform.filter(id=transform.id).first()
            if transform_exists is None:
                transform.save()
                self._logging_message += f"created Transform('{transform.uid[:8]}')"
                transform_exists = transform
            else:
                self._logging_message += f"loaded Transform('{transform.uid[:8]}')"
            self._transform = transform_exists

        if new_run is None:  # for notebooks, default to loading latest runs
            new_run = False if self._transform.type == "notebook" else True  # type: ignore

        run = None
        if not new_run:  # try loading latest run by same user
            run = (
                Run.filter(
                    transform=self._transform, created_by_id=ln_setup.settings.user.id
                )
                .order_by("-created_at")
                .first()
            )
            if run is not None:  # loaded latest run
                run.started_at = datetime.now(timezone.utc)  # update run time
                self._logging_message += f", started Run('{run.uid[:8]}') at {format_field_value(run.started_at)}"

        if run is None:  # create new run
            run = Run(
                transform=self._transform,
                params=params,
            )
            run.started_at = datetime.now(timezone.utc)
            self._logging_message += f", started new Run('{run.uid[:8]}') at {format_field_value(run.started_at)}"
        # can only determine at ln.finish() if run was consecutive in
        # interactive session, otherwise, is consecutive
        run.is_consecutive = True if is_run_from_ipython else None
        # need to save in all cases
        run.save()
        if params is not None:
            run.params.add_values(params)
            self._logging_message += "\n→ params: " + " ".join(
                f"{key}='{value}'" for key, value in params.items()
            )
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
        key = Path(path).name
        if isinstance(path, (Path, PurePath)):
            path_str = path.as_posix()  # type: ignore
        else:
            path_str = str(path)
        if path_str.endswith("Untitled.ipynb"):
            raise RuntimeError("Please rename your notebook before tracking it")
        if path_str.startswith("/fileId="):
            name = get_notebook_name_colab()
            key = f"{name}.ipynb"
        else:
            import nbproject

            try:
                nbproject_title = nbproject.meta.live.title
            except IndexError:
                raise NotebookNotSaved(
                    "The notebook is not saved, please save the notebook and"
                    " rerun ``"
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
        def get_key_clashing_message(transform: Transform, key: str) -> str:
            update_key_note = message_update_key_in_version_family(
                suid=transform.stem_uid,
                existing_key=transform.key,
                new_key=key,
                registry="Transform",
            )
            return (
                f'Filename "{key}" clashes with the existing key "{transform.key}" for uid "{transform.uid[:-4]}...."\n\nEither init a new transform with a new uid:\n\n'
                f'ln.track("{ids.base62_12()}0000)"\n\n{update_key_note}'
            )

        # make a new transform record
        if transform is None:
            if uid is None:
                uid = f"{stem_uid}{get_uid_ext(version)}"
            # let's query revises so that we can pass it to the constructor and use it for error handling
            revises = (
                Transform.filter(uid__startswith=uid[:-4], is_latest=True)
                .order_by("-created_at")
                .first()
            )
            # note that here we're not passing revises because we're not querying it
            # hence, we need to do a revision family lookup based on key
            # hence, we need key to be not None
            assert key is not None  # noqa: S101
            raise_update_context = False
            try:
                transform = Transform(
                    uid=uid,
                    version=version,
                    name=name,
                    key=key,
                    reference=transform_ref,
                    reference_type=transform_ref_type,
                    type=transform_type,
                    revises=revises,
                ).save()
            except InconsistentKey:
                raise_update_context = True
            if raise_update_context:
                raise UpdateContext(get_key_clashing_message(revises, key))
            self._logging_message += f"created Transform('{transform.uid[:8]}')"
        else:
            uid = transform.uid
            # transform was already saved via `finish()`
            transform_was_saved = (
                transform._source_code_artifact_id is not None
                or transform.source_code is not None
            )
            # check whether the transform.key is consistent
            if transform.key != key:
                raise UpdateContext(get_key_clashing_message(transform, key))
            elif transform.name != name:
                transform.name = name
                transform.save()
                self._logging_message += (
                    "updated transform name, "  # white space on purpose
                )
            elif (
                transform.created_by_id != ln_setup.settings.user.id
                and not transform_was_saved
            ):
                raise UpdateContext(
                    f'{transform.created_by.name} ({transform.created_by.handle}) already works on this draft {transform.type}.\n\nPlease create a revision via `ln.track("{uid[:-4]}{increment_base62(uid[-4:])}")` or a new transform with a *different* filename and `ln.track("{ids.base62_12()}0000")`.'
                )
            # check whether transform source code was already saved
            if transform_was_saved:
                bump_revision = False
                if is_run_from_ipython:
                    bump_revision = True
                else:
                    hash, _ = hash_file(self._path)  # ignore hash_type for now
                    if transform.hash is not None:
                        condition = hash != transform.hash
                    else:
                        condition = hash != transform._source_code_artifact.hash
                    if condition:
                        bump_revision = True
                    else:
                        self._logging_message += (
                            f"loaded Transform('{transform.uid[:8]}')"
                        )
                if bump_revision:
                    change_type = (
                        "Re-running saved notebook"
                        if is_run_from_ipython
                        else "Source code changed"
                    )
                    raise UpdateContext(
                        f"{change_type}, bump revision by setting:\n\n"
                        f'ln.track("{uid[:-4]}{increment_base62(uid[-4:])}")'
                    )
            else:
                self._logging_message += f"loaded Transform('{transform.uid[:8]}')"
        self._transform = transform

    def finish(self, ignore_non_consecutive: None | bool = None) -> None:
        """Finish a tracked run.

        - writes a timestamp: `run.finished_at`
        - saves the source code: `transform.source_code`

        When called in the last cell of a notebook:

        - prompts for user input if not consecutively executed
        - requires to save the notebook in your editor right before
        - saves a run report: `run.report`

        Args:
            ignore_non_consecutive: Whether to ignore if a notebook was non-consecutively executed.

        Examples:

            >>> import lamindb as ln
            >>> ln.track()
            >>> # do things while tracking data lineage
            >>> ln.finish()

        See Also:
            `lamin save script.py` or `lamin save notebook.ipynb` → `docs </cli#lamin-save>`__

        """
        from lamindb._finish import save_context_core

        def get_seconds_since_modified(filepath) -> float:
            return datetime.now().timestamp() - filepath.stat().st_mtime

        def get_shortcut() -> str:
            import platform

            return "CMD + s" if platform.system() == "Darwin" else "CTRL + s"

        if self.run is None:
            raise TrackNotCalled("Please run `ln.track()` before `ln.finish()`")
        if self._path is None:
            if self.run.transform.type in {"script", "notebook"}:
                raise ValueError(
                    "Transform type is not allowed to be 'script' or 'notebook' because `context._path` is `None`."
                )
            self.run.finished_at = datetime.now(timezone.utc)
            self.run.save()
            # nothing else to do
            return None
        if is_run_from_ipython:  # notebooks
            import nbproject

            # it might be that the user modifies the title just before ln.finish()
            if (nbproject_title := nbproject.meta.live.title) != self.transform.name:
                self.transform.name = nbproject_title
                self.transform.save()
            if get_seconds_since_modified(self._path) > 2 and not ln_setup._TESTING:
                raise NotebookNotSaved(
                    f"Please save the notebook in your editor (shortcut `{get_shortcut()}`) right before calling `ln.finish()`"
                )
        save_context_core(
            run=self.run,
            transform=self.run.transform,
            filepath=self._path,
            finished_at=True,
            ignore_non_consecutive=ignore_non_consecutive,
        )


context = Context()
