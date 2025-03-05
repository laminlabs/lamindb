from __future__ import annotations

import builtins
import hashlib
import signal
import sys
import threading
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
from django.db.models import Func, IntegerField, Q
from lamin_utils import logger
from lamindb_setup.core import deprecated
from lamindb_setup.core.hashing import hash_file

from lamindb.base import ids
from lamindb.base.ids import base62_12
from lamindb.models import Run, Transform, format_field_value

from ..core._settings import settings
from ..errors import (
    InconsistentKey,
    InvalidArgument,
    TrackNotCalled,
    UpdateContext,
)
from ..models._is_versioned import bump_version as bump_version_function
from ..models._is_versioned import (
    increment_base62,
    message_update_key_in_version_family,
)
from ._sync_git import get_transform_reference_from_git_repo
from ._track_environment import track_environment

if TYPE_CHECKING:
    from lamindb_setup.core.types import UPathStr

    from lamindb.base.types import TransformType
    from lamindb.models import Project

is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)

msg_path_failed = "failed to infer notebook path.\nfix: pass `path` to `ln.track()`"


def get_uid_ext(version: str) -> str:
    from lamin_utils._base62 import encodebytes

    # merely zero-padding the nbproject version such that the base62 encoding is
    # at least 4 characters long doesn't yields sufficiently diverse hashes and
    # leads to collisions; it'd be nice because the uid_ext would be ordered
    return encodebytes(hashlib.md5(version.encode()).digest())[:4]  # noqa: S324


def get_notebook_path() -> Path:
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
    return Path(path)


# from https://stackoverflow.com/questions/61901628
def get_notebook_key_colab() -> str:
    from socket import gethostbyname, gethostname  # type: ignore

    from requests import get  # type: ignore

    ip = gethostbyname(gethostname())  # 172.28.0.12
    try:
        key = get(f"http://{ip}:9000/api/sessions").json()[0]["name"]  # noqa: S113
        key = f"colab/{key}"
    except Exception:
        logger.warning(
            "could not get notebook key from Google Colab, using: colab/notebook.ipynb"
        )
        key = "colab/notebook.ipynb"
    return key


def pretty_pypackages(dependencies: dict) -> str:
    deps_list = []
    for pkg, ver in dependencies.items():
        if ver != "":
            deps_list.append(pkg + f"=={ver}")
        else:
            deps_list.append(pkg)
    deps_list.sort()
    return " ".join(deps_list)


class LogStreamHandler:
    def __init__(self, log_stream, file):
        self.log_stream = log_stream
        self.file = file

    def write(self, data):
        self.log_stream.write(data)
        self.file.write(data)
        self.file.flush()

    def flush(self):
        self.log_stream.flush()
        self.file.flush()


class LogStreamTracker:
    def __init__(self):
        self.original_stdout = None
        self.original_stderr = None
        self.log_file = None
        self.original_excepthook = sys.excepthook
        self.is_cleaning_up = False

    def start(self, run: Run):
        self.original_stdout = sys.stdout
        self.original_stderr = sys.stderr
        self.run = run
        self.log_file_path = (
            ln_setup.settings.cache_dir / f"run_logs_{self.run.uid}.txt"
        )
        self.log_file = open(self.log_file_path, "w")
        sys.stdout = LogStreamHandler(self.original_stdout, self.log_file)
        sys.stderr = LogStreamHandler(self.original_stderr, self.log_file)
        # handle signals
        # signal should be used only in the main thread, otherwise
        # ValueError: signal only works in main thread of the main interpreter
        if threading.current_thread() == threading.main_thread():
            signal.signal(signal.SIGTERM, self.cleanup)
            signal.signal(signal.SIGINT, self.cleanup)
        # handle exceptions
        sys.excepthook = self.handle_exception

    def finish(self):
        if self.original_stdout:
            sys.stdout = self.original_stdout
            sys.stderr = self.original_stderr
            self.log_file.close()

    def cleanup(self, signo=None, frame=None):
        from lamindb._finish import save_run_logs

        if self.original_stdout and not self.is_cleaning_up:
            self.is_cleaning_up = True
            if signo is not None:
                signal_msg = f"\nProcess terminated by signal {signo} ({signal.Signals(signo).name})\n"
                if frame:
                    signal_msg += (
                        f"Frame info:\n{''.join(traceback.format_stack(frame))}"
                    )
                self.log_file.write(signal_msg)
            sys.stdout = self.original_stdout
            sys.stderr = self.original_stderr
            self.log_file.flush()
            self.log_file.close()
            save_run_logs(self.run, save_run=True)

    def handle_exception(self, exc_type, exc_value, exc_traceback):
        if not self.is_cleaning_up:
            error_msg = f"{''.join(traceback.format_exception(exc_type, exc_value, exc_traceback))}"
            if self.log_file.closed:
                self.log_file = open(self.log_file_path, "a")
            self.log_file.write(error_msg)
            self.log_file.flush()
            self.cleanup()
        self.original_excepthook(exc_type, exc_value, exc_traceback)


class Context:
    """Run context.

    Enables convenient data lineage tracking by managing a transform & run
    upon :meth:`~lamindb.core.Context.track` & :meth:`~lamindb.core.Context.finish`.

    Guide: :doc:`/track`

    Examples:

        Is typically used via the global :class:`~lamindb.context` object via `ln.track()` and `ln.finish()`:

        >>> import lamindb as ln
        >>> ln.track()
        >>> # do things
        >>> ln.finish()

    """

    def __init__(self):
        self._uid: str | None = None
        self._description: str | None = None
        self._version: str | None = None
        self._transform: Transform | None = None
        self._run: Run | None = None
        self._path: Path | None = None
        """A local path to the script that's running."""
        self._project: Project | None = None
        self._logging_message_track: str = ""
        self._logging_message_imports: str = ""
        self._stream_tracker: LogStreamTracker = LogStreamTracker()
        self._is_finish_retry: bool = False

    @property
    def transform(self) -> Transform | None:
        """Managed transform of context."""
        return self._transform

    @property
    def description(self) -> str | None:
        """`description` argument for `context.transform`."""
        return self._description

    @description.setter
    def description(self, value: str | None):
        self._description = value

    @property
    @deprecated(new_name="description")
    def name(self) -> str | None:
        return self._description

    @name.setter
    def name(self, value: str | None):
        self._description = value

    @property
    def uid(self) -> str | None:
        """`uid` argument for `context.transform`."""
        return self._uid

    @uid.setter
    def uid(self, value: str | None):
        self._uid = value

    @property
    def version(self) -> str | None:
        """`version` argument for `context.transform`."""
        return self._version

    @version.setter
    def version(self, value: str | None):
        self._version = value

    @property
    def project(self) -> Project | None:
        """Project to label entities created during the run."""
        return self._project

    @property
    def run(self) -> Run | None:
        """Managed run of context."""
        return self._run

    def track(
        self,
        transform: str | Transform | None = None,
        *,
        project: str | None = None,
        params: dict | None = None,
        new_run: bool | None = None,
        path: str | None = None,
    ) -> None:
        """Track a global run of your Python session.

        - sets :attr:`~lamindb.core.Context.transform` &
          :attr:`~lamindb.core.Context.run` by creating or loading `Transform` &
          `Run` records
        - saves Python environment as a `requirements.txt` file: `run.environment`

        If :attr:`~lamindb.core.Settings.sync_git_repo` is set, checks whether a
        script-like transform exists in a git repository and links it.

        Args:
            transform: A transform `uid` or record. If `None`, creates a `uid`.
            project: A project `name` or `uid` for labeling entities created during the run.
            params: A dictionary of parameters to track for the run.
            new_run: If `False`, loads the latest run of transform
                (default notebook), if `True`, creates new run (default non-notebook).
            path: Filepath of notebook or script. Only needed if it can't be
                automatically detected.

        Examples:

            To track the run of a notebook or script, call:

            >>> ln.track()

            If you want to ensure a single version history across renames of the notebook or script, pass the auto-generated `uid` that you'll find in the logs:

            >>> ln.track("Onv04I53OgtT0000")  # example uid, the last four characters encode the version of the transform

        """
        from lamindb.models import Project

        if project is not None:
            project_record = Project.filter(
                Q(name=project) | Q(uid=project)
            ).one_or_none()
            if project_record is None:
                raise InvalidArgument(
                    f"Project '{project}' not found, either create it with `ln.Project(name='...').save()` or fix typos."
                )
            self._project = project_record
        self._logging_message_track = ""
        self._logging_message_imports = ""
        if transform is not None and isinstance(transform, str):
            self.uid = transform
            transform = None
        self._path = None
        if transform is None:
            description = None
            if is_run_from_ipython:
                self._path, description = self._track_notebook(path_str=path)
                transform_type = "notebook"
                transform_ref = None
                transform_ref_type = None
            else:
                (
                    self._path,
                    transform_type,
                    transform_ref,
                    transform_ref_type,
                ) = self._track_source_code(path=path)
            if description is None:
                description = self._description
            # temporarily until the hub displays the key by default
            # populate the description with the filename again
            if description is None:
                description = self._path.name
            self._create_or_load_transform(
                description=description,
                transform_ref=transform_ref,
                transform_ref_type=transform_ref_type,
                transform_type=transform_type,  # type: ignore
            )
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
                self._logging_message_track += f"created Transform('{transform.uid}')"
                transform_exists = transform
            else:
                self._logging_message_track += f"loaded Transform('{transform.uid}')"
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
                self._logging_message_track += f", re-started Run('{run.uid[:8]}...') at {format_field_value(run.started_at)}"

        if run is None:  # create new run
            run = Run(  # type: ignore
                transform=self._transform,
                params=params,
            )
            run.started_at = datetime.now(timezone.utc)
            self._logging_message_track += f", started new Run('{run.uid[:8]}...') at {format_field_value(run.started_at)}"
        # can only determine at ln.finish() if run was consecutive in
        # interactive session, otherwise, is consecutive
        run.is_consecutive = True if is_run_from_ipython else None
        # need to save in all cases
        run.save()
        if params is not None:
            run.params.add_values(params)
            self._logging_message_track += "\n→ params: " + ", ".join(
                f"{key}={value}" for key, value in params.items()
            )
        self._run = run
        track_environment(run)
        if self.project is not None:
            # to update a potential project link
            # is only necessary if transform is loaded rather than newly created
            # can be optimized by checking whether the transform is loaded, but it typically is
            transform.save()
        log_to_file = None
        if log_to_file is None:
            log_to_file = self.transform.type != "notebook"
        if log_to_file:
            self._stream_tracker.start(run)
        logger.important(self._logging_message_track)
        if self._logging_message_imports:
            logger.important(self._logging_message_imports)

    def _track_source_code(
        self,
        *,
        path: UPathStr | None,
    ) -> tuple[Path, str, str, str]:
        # for `.py` files, classified as "script"
        # for `.Rmd` and `.qmd` files, which we classify
        # as "notebook" because they typically come with an .html run report
        if path is None:
            import inspect

            frame = inspect.stack()[2]
            module = inspect.getmodule(frame[0])
            # None for interactive session
            if module is None:
                raise NotImplementedError(
                    "Interactive sessions are not yet supported to be tracked."
                )
            path = Path(module.__file__)
        else:
            path = Path(path)
        # for Rmd and qmd, we could also extract the title
        # we don't do this for now as we're setting the title upon `ln.finish()` or `lamin save`
        # by extracting it from the html while cleaning it: see clean_r_notebook_html()
        # also see the script_to_notebook() in the CLI _load.py where the title is extracted
        # from the source code YAML and updated with the transform description
        transform_type = "notebook" if path.suffix in {".Rmd", ".qmd"} else "script"
        reference = None
        reference_type = None
        if settings.sync_git_repo is not None:
            reference = get_transform_reference_from_git_repo(path)
            reference_type = "url"
        return path, transform_type, reference, reference_type

    def _track_notebook(
        self,
        *,
        path_str: str | None,
    ) -> tuple[Path, str | None]:
        if path_str is None:
            path = get_notebook_path()
        else:
            path = Path(path_str)
        description = None
        path_str = path.as_posix()
        if path_str.endswith("Untitled.ipynb"):
            raise RuntimeError("Please rename your notebook before tracking it")
        if path_str.startswith("/fileId="):
            logger.warning("tracking on Google Colab is experimental")
            path_str = get_notebook_key_colab()
            path = Path(path_str)
        else:
            import nbproject

            try:
                nbproject_title = nbproject.meta.live.title
            except IndexError:
                # notebook is not saved
                pass
            if nbproject_title is not None:
                description = nbproject_title
            # log imported python packages
            try:
                from nbproject.dev._pypackage import infer_pypackages

                nb = nbproject.dev.read_notebook(path_str)
                self._logging_message_imports += (
                    "notebook imports:"
                    f" {pretty_pypackages(infer_pypackages(nb, pin_versions=True))}"
                )
            except Exception:
                logger.debug("inferring imported packages failed")
                pass
        return path, description

    def _create_or_load_transform(
        self,
        *,
        description: str,
        transform_ref: str | None = None,
        transform_ref_type: str | None = None,
        transform_type: TransformType = None,
    ):
        def get_key_clashing_message(transform: Transform, key: str) -> str:
            update_key_note = message_update_key_in_version_family(
                suid=transform.stem_uid,
                existing_key=transform.key,
                new_key=key,
                registry="Transform",
            )
            return (
                f'Filepath "{key}" clashes with the existing key "{transform.key}" for uid "{transform.uid[:-4]}...."\n\nEither init a new transform with a new uid:\n\n'
                f'ln.track("{ids.base62_12()}0000")\n\n{update_key_note}'
            )

        revises = None
        # the user did not pass the uid
        if self.uid is None:

            class SlashCount(Func):
                template = "LENGTH(%(expressions)s) - LENGTH(REPLACE(%(expressions)s, '/', ''))"
                output_field = IntegerField()

            # we need to traverse from greater depth to shorter depth so that we match better matches first
            transforms = (
                Transform.filter(key__endswith=self._path.name, is_latest=True)
                .annotate(slash_count=SlashCount("key"))
                .order_by("-slash_count")
            )
            uid = f"{base62_12()}0000"
            key = self._path.name
            target_transform = None
            hash, _ = hash_file(self._path)
            if len(transforms) != 0:
                message = ""
                found_key = False
                for aux_transform in transforms:
                    if aux_transform.key in self._path.as_posix():
                        key = aux_transform.key
                        if (
                            # has to be the same user
                            aux_transform.created_by_id == ln_setup.settings.user.id
                            and (
                                # if the transform source code wasn't yet saved
                                aux_transform.source_code is None
                                # if the transform source code is unchanged
                                # if aux_transform.type == "notebook", we anticipate the user makes changes to the notebook source code
                                # in an interactive session, hence we *pro-actively bump* the version number by setting `revises`
                                # in the second part of the if condition even though the source code is unchanged at point of running track()
                                or (
                                    aux_transform.hash == hash
                                    and aux_transform.type != "notebook"
                                )
                            )
                        ):
                            uid = aux_transform.uid
                            target_transform = aux_transform
                        else:
                            uid = f"{aux_transform.uid[:-4]}{increment_base62(aux_transform.uid[-4:])}"
                            message = f"there already is a {aux_transform.type} with `key` '{aux_transform.key}'"
                            if (
                                aux_transform.hash == hash
                                and aux_transform.type == "notebook"
                            ):
                                message += " -- anticipating changes"
                            elif aux_transform.hash != hash:
                                message += ""  # could log "source code changed", but this seems too much
                            elif (
                                aux_transform.created_by_id != ln_setup.settings.user.id
                            ):
                                message += f" -- {aux_transform.created_by.handle} already works on this draft"
                            message += f", creating new version '{uid}'"
                            revises = aux_transform
                        found_key = True
                        break
                if not found_key:
                    plural_s = "s" if len(transforms) > 1 else ""
                    transforms_str = "\n".join(
                        [
                            f"    {transform.uid} → {transform.key}"
                            for transform in transforms
                        ]
                    )
                    message = f"ignoring transform{plural_s} with same filename:\n{transforms_str}"
                if message != "":
                    logger.important(message)
            self.uid, transform = uid, target_transform
        # the user did pass the uid
        else:
            transform = Transform.filter(uid=self.uid).one_or_none()
            if transform is not None:
                if transform.key not in self._path.as_posix():
                    n_parts = len(Path(transform.key).parts)
                    last_path_elements = (
                        Path(*self._path.parts[-n_parts:]).as_posix()
                        if n_parts > 0
                        else ""
                    )
                    raise UpdateContext(
                        get_key_clashing_message(transform, last_path_elements)
                    )
                key = transform.key  # type: ignore
            else:
                key = self._path.name
        if self.version is not None:
            # test inconsistent version passed
            if (
                transform is not None
                and transform.version is not None  # type: ignore
                and self.version != transform.version  # type: ignore
            ):
                raise SystemExit(
                    f"✗ please pass consistent version: ln.context.version = '{transform.version}'"  # type: ignore
                )
            # test whether version was already used for another member of the family
            suid, vuid = (self.uid[:-4], self.uid[-4:])
            transform = Transform.filter(
                uid__startswith=suid, version=self.version
            ).one_or_none()
            if transform is not None and vuid != transform.uid[-4:]:
                better_version = bump_version_function(self.version)
                raise SystemExit(
                    f"✗ version '{self.version}' is already taken by Transform('{transform.uid}'); please set another version, e.g., ln.context.version = '{better_version}'"
                )
        # make a new transform record
        if transform is None:
            assert key is not None  # noqa: S101
            raise_update_context = False
            try:
                transform = Transform(  # type: ignore
                    uid=self.uid,
                    version=self.version,
                    description=description,
                    key=key,
                    reference=transform_ref,
                    reference_type=transform_ref_type,
                    type=transform_type,
                ).save()
            except InconsistentKey:
                raise_update_context = True
            if raise_update_context:
                if revises is None:
                    revises = (
                        Transform.filter(uid__startswith=self.uid[:-4], is_latest=True)
                        .order_by("-created_at")
                        .first()
                    )
                raise UpdateContext(get_key_clashing_message(revises, key))
            self._logging_message_track += f"created Transform('{transform.uid}')"
        else:
            uid = transform.uid
            # transform was already saved via `finish()`
            transform_was_saved = transform.source_code is not None
            # check whether the transform.key is consistent
            if transform.key != key:
                raise UpdateContext(get_key_clashing_message(transform, key))
            elif transform.description != description and description is not None:
                transform.description = description
                transform.save()
                self._logging_message_track += (
                    "updated transform description, "  # white space on purpose
                )
            elif (
                transform.created_by_id != ln_setup.settings.user.id
                and not transform_was_saved
            ):
                raise UpdateContext(
                    f'{transform.created_by.name} ({transform.created_by.handle}) already works on this draft {transform.type}.\n\nPlease create a revision via `ln.track("{uid[:-4]}{increment_base62(uid[-4:])}")` or a new transform with a *different* filedescription and `ln.track("{ids.base62_12()}0000")`.'
                )
            # check whether transform source code was already saved
            if transform_was_saved:
                bump_revision = False
                if transform.type == "notebook":
                    # we anticipate the user makes changes to the notebook source code
                    # in an interactive session, hence we pro-actively bump the version number
                    bump_revision = True
                else:
                    hash, _ = hash_file(self._path)  # ignore hash_type for now
                    if hash != transform.hash:
                        bump_revision = True
                    else:
                        self._logging_message_track += (
                            f"loaded Transform('{transform.uid}')"
                        )
                if bump_revision:
                    change_type = (
                        "re-running notebook with already-saved source code"
                        if transform.type == "notebook"
                        else "source code changed"
                    )
                    raise UpdateContext(
                        f'✗ {change_type}, please update the `uid` argument in `track()` to "{uid[:-4]}{increment_base62(uid[-4:])}"'
                    )
            else:
                self._logging_message_track += f"loaded Transform('{transform.uid}')"
        self._transform = transform

    def finish(self, ignore_non_consecutive: None | bool = None) -> None:
        """Finish a tracked run.

        - writes a timestamp: `run.finished_at`
        - saves the source code: `transform.source_code`
        - saves a run report: `run.report`

        When called in the last cell of a notebook:

        - prompts to save the notebook in your editor right before
        - prompts for user input if not consecutively executed

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
        from lamindb._finish import (
            save_context_core,
        )

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
        return_code = save_context_core(
            run=self.run,
            transform=self.run.transform,
            filepath=self._path,
            finished_at=True,
            ignore_non_consecutive=ignore_non_consecutive,
            is_retry=self._is_finish_retry,
        )
        if return_code == "retry":
            self._is_finish_retry = True
            return None
        if self.transform.type != "notebook":
            self._stream_tracker.finish()
        # reset the context attributes so that somebody who runs `track()` after finish
        # starts fresh
        self._uid = None
        self._run = None
        self._transform = None
        self._version = None
        self._description = None


context = Context()
