from __future__ import annotations

import builtins
import hashlib
import os
import signal
import sys
import threading
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, TextIO

import lamindb_setup as ln_setup
from django.db.models import Q
from lamin_utils._logger import logger
from lamindb_setup.core.django import _is_running_in_marimo

from .._secret_redaction import (
    REDACTED_SECRET_VALUE,
    is_sensitive_param_key,
    is_sensitive_param_value,
)
from ..errors import InvalidArgument, TrackNotCalled
from ..models import Run, SQLRecord, Transform, format_field_value
from ..models._feature_manager import infer_convert_dtype_key_value
from ._settings import settings
from ._sync_git import get_transform_reference_from_git_repo
from ._track_environment import track_python_environment

if TYPE_CHECKING:
    from types import FrameType, TracebackType

    from lamindb.base.types import TransformKind
    from lamindb.models import Artifact, Branch, Project, Space


is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)

msg_path_failed = "failed to infer notebook path.\nfix: pass `path` to `ln.track()`"


def get_key_from_module(caller_module: str) -> str:
    if "." in caller_module:
        key_from_module = f"pypackages/{caller_module.replace('.', '/')}.py"
    else:
        key_from_module = None
    return key_from_module


def detect_and_process_source_code_file(
    *,
    path: str | Path | None,
    transform_kind: TransformKind | None = None,
    infer_reference: bool = True,
) -> tuple[Path, TransformKind, str, str, str | None]:
    """Track source code file and determine transform metadata.

    For `.py` files, classified as "script".
    For `.Rmd` and `.qmd` files, classified as "notebook" because they
    typically come with an .html run report.

    Package vs script criterion: source code is part of a **package** if the
    caller's module name contains at least one `.` (module nesting goes beyond
    the filename). Otherwise it is a **script** (module nesting stops at the
    filename, e.g. `__main__`, `__mp_main__`, or a single top-level name).

    Args:
        path: Path to the source code file. If None, infers from call stack.

    Returns:
        Tuple of (path, transform_kind, reference, reference_type, key_from_module).
        - path: Path object to the source file
        - transform_kind: "script" or "notebook"
        - reference: Git reference URL if sync_git_repo is set, else None
        - reference_type: "url" if reference exists, else None
        - key_from_module: If caller is part of a package (`.` in __name__),
          `pypackages/module/path/to/file.py`; else None (key will be computed from dev_dir or path.name).

    Raises:
        NotImplementedError: If path cannot be determined from call stack.
    """
    # for `.py` files, classified as "script"
    # for `.Rmd` and `.qmd` files, which we classify
    # as "notebook" because they typically come with an .html run report
    key_from_module: str | None = None
    if path is None:
        import inspect

        frame = inspect.stack()[2]
        path_str = frame[1]
        if not path_str or path_str.startswith("<"):
            raise NotImplementedError(
                "Cannot determine valid file path, pass manually via path (interactive sessions not yet supported)"
            )
        path = Path(path_str)
        # package vs script: nesting beyond filename makes the file part of a python package
        caller_module = frame[0].f_globals.get("__name__", "__main__")
        key_from_module = get_key_from_module(caller_module)
    else:
        path = Path(path)
    # for Rmd and qmd, we could also extract the title
    # we don't do this for now as we're setting the title upon `ln.finish()` or `lamin save`
    # by extracting it from the html while cleaning it: see clean_r_notebook_html()
    # also see the script_to_notebook() in the CLI _load.py where the title is extracted
    # from the source code YAML and updated with the transform description
    # note that ipynb notebooks are handled in a separate function (_track_notebook())
    if transform_kind is None:
        transform_kind = "notebook" if path.suffix in {".Rmd", ".qmd"} else "script"
    reference = None
    reference_type = None
    if (
        infer_reference
        and settings.sync_git_repo is not None
        and path.suffix != ".ipynb"
    ):
        reference = get_transform_reference_from_git_repo(path)
        reference_type = "url"
    return path, transform_kind, reference, reference_type, key_from_module


def get_uid_ext(version: str) -> str:
    from lamin_utils._base62 import encodebytes

    # merely zero-padding the nbproject version such that the base62 encoding is
    # at least 4 characters long doesn't yields sufficiently diverse hashes and
    # leads to collisions; it'd be nice because the uid_ext would be ordered
    return encodebytes(hashlib.md5(version.encode()).digest())[:4]  # noqa: S324


def get_marimo_notebook_path() -> str | None:
    if not _is_running_in_marimo():
        return None

    from marimo._runtime import context as marimo_runtime_context

    if not hasattr(marimo_runtime_context, "safe_get_context"):
        raise RuntimeError(
            "marimo version is incompatible with lamindb: "
            "marimo._runtime.context lacks safe_get_context"
        )

    ctx = marimo_runtime_context.safe_get_context()
    if ctx is None:
        return None
    df_mode = ctx.marimo_config.get("display", {}).get("dataframes")
    if df_mode != "plain":
        logger.warning(
            "marimo's `Dataframe viewer` will make tables in your "
            "run report appear as truncated pngs - fix: "
            "in the top left corner of the UI, go to Settings → User settings → Packages & Data → Data and set `Dataframe viewer` to `plain`"
        )
    return getattr(ctx, "filename", None)


def get_notebook_path() -> tuple[Path, str]:
    marimo_path = get_marimo_notebook_path()
    if marimo_path is not None:
        return Path(marimo_path), "marimo"

    from nbproject.dev._jupyter_communicate import (
        notebook_path as get_notebook_path,
    )

    path = None
    try:
        path, env = get_notebook_path(return_env=True)
    except ValueError as ve:
        raise ve
    except Exception as error:
        raise RuntimeError(msg_path_failed) from error
    if path is None:
        raise RuntimeError(msg_path_failed) from None
    return Path(path), env


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


def get_cli_call() -> tuple[str, str] | None:
    """Returns (tool_name, args) when invoked as a script with CLI arguments.

    Returns None if not run as a script (e.g., in Jupyter, interactive shell)
    or when no arguments were passed.
    """
    if len(sys.argv) > 1 and sys.argv[0] and not is_run_from_ipython:
        return Path(sys.argv[0]).name, " ".join(sys.argv[1:])
    return None


def pretty_pypackages(dependencies: dict) -> str:
    deps_list = []
    for pkg, ver in dependencies.items():
        if ver != "":
            deps_list.append(pkg + f"=={ver}")
        else:
            deps_list.append(pkg)
    deps_list.sort()
    return " ".join(deps_list)


def last_non_empty_r_block(line: str) -> str:
    for block in reversed(line.split("\r")):
        if block:
            return block
    return ""


class LogStreamHandler:
    def __init__(self, log_stream: TextIO, file: TextIO, use_buffer: bool):
        self.log_stream = log_stream
        self.file = file

        self._buffer = ""
        self._use_buffer = use_buffer

    def write(self, data: str) -> int:
        data_length = len(data)

        self.log_stream.write(data)
        if self.file.closed:
            return data_length

        if not self._use_buffer:
            self.file.write(data)
            self.file.flush()
            return data_length

        self._buffer += data
        # write only the last part of a line with carriage returns
        while "\n" in self._buffer:
            if self.file.closed:
                return data_length
            line, self._buffer = self._buffer.split("\n", 1)
            self.file.write(last_non_empty_r_block(line) + "\n")
            self.file.flush()

        return data_length

    def flush(self):
        self.log_stream.flush()
        if not self.file.closed:
            self.file.flush()

    # https://laminlabs.slack.com/archives/C07DB677JF6/p1759423901926139
    # other tracking frameworks like W&B use our output stream and expect
    # certain functions like isatty to be available
    def isatty(self) -> bool:
        return False

    # .flush is sometimes (in jupyter etc.) called after every .write
    # this needs to be called only at the end
    def flush_buffer(self):
        if not self.file.closed and self._buffer:
            self.file.write(last_non_empty_r_block(self._buffer))
            self._buffer = ""
        self.flush()


class LogStreamTracker:
    def __init__(self):
        self.original_stdout = None
        self.original_stderr = None
        self.log_file = None
        self.is_cleaning_up = False
        self.original_excepthook: Callable[
            [type[BaseException], BaseException, TracebackType | None], Any
        ] = sys.excepthook

        self.original_signal_handlers: dict[
            signal.Signals, Callable[[int, FrameType | None], Any] | int
        ] = {}
        if threading.current_thread() == threading.main_thread():
            self.original_signal_handlers[signal.SIGTERM] = signal.getsignal(
                signal.SIGTERM
            )
            self.original_signal_handlers[signal.SIGINT] = signal.getsignal(
                signal.SIGINT
            )

    def start(self, run: Run):
        self.original_stdout = sys.stdout
        self.original_stderr = sys.stderr
        self.run = run
        self.log_file_path = (
            ln_setup.settings.cache_dir / f"run_logs_{self.run.uid}.txt"
        )
        self.log_file = open(self.log_file_path, "w", encoding="utf-8")
        # the instance that's connected is important information
        self.log_file.write(
            f"\x1b[92m→\x1b[0m connected lamindb: {ln_setup.settings.instance.slug}\n"
        )
        # use buffering for correct handling of carriage returns
        sys.stdout = LogStreamHandler(
            self.original_stdout, self.log_file, use_buffer=True
        )
        # write evrything immediately in stderr
        sys.stderr = LogStreamHandler(
            self.original_stderr, self.log_file, use_buffer=False
        )
        # handle signals
        # signal should be used only in the main thread, otherwise
        # ValueError: signal only works in main thread of the main interpreter
        if threading.current_thread() == threading.main_thread():
            signal.signal(signal.SIGTERM, self.cleanup)
            signal.signal(signal.SIGINT, self.cleanup)
        # handle exceptions
        sys.excepthook = self.handle_exception
        # reset handler for lamin logger because sys.stdout has been replaced
        logger.set_handler()

    def finish(self):
        if self.original_stdout:
            getattr(sys.stdout, "flush_buffer", sys.stdout.flush)()
            sys.stderr.flush()
            sys.stdout = self.original_stdout
            sys.stderr = self.original_stderr
            if not self.log_file.closed:
                self.log_file.close()
            # reset handler for lamin logger because sys.stdout has been replaced
            logger.set_handler()

    def cleanup(self, signo=None, frame=None):
        try:
            from .._finish import save_run_logs

            if self.original_stdout and not self.is_cleaning_up:
                self.is_cleaning_up = True
                if signo is not None:
                    if self.log_file.closed:
                        self.log_file = open(self.log_file_path, "a", encoding="utf-8")
                    getattr(sys.stdout, "flush_buffer", sys.stdout.flush)()
                    sys.stderr.flush()
                    signal_msg = f"\nProcess terminated by signal {signo} ({signal.Signals(signo).name})\n"
                    if frame:
                        signal_msg += (
                            f"Frame info:\n{''.join(traceback.format_stack(frame))}"
                        )
                    self.log_file.write(signal_msg)
                    self.log_file.flush()
                    self.run._status_code = 2  # aborted
                else:
                    self.run._status_code = 1  # errored
                self.run.finished_at = datetime.now(timezone.utc)
                sys.stdout = self.original_stdout
                sys.stderr = self.original_stderr
                if not self.log_file.closed:
                    self.log_file.close()
                save_run_logs(self.run, save_run=True)
                # reset handler for lamin logger because sys.stdout has been replaced
                logger.set_handler()
        except:  # noqa: E722, S110
            pass
        finally:
            if signo is not None and signo in self.original_signal_handlers:
                original_handler = self.original_signal_handlers[signo]
                if callable(original_handler):
                    original_handler(signo, frame)

    def handle_exception(self, exc_type, exc_value, exc_traceback):
        try:
            if self.original_stdout and not self.is_cleaning_up:
                if self.log_file.closed:
                    self.log_file = open(self.log_file_path, "a", encoding="utf-8")
                getattr(sys.stdout, "flush_buffer", sys.stdout.flush)()
                sys.stderr.flush()
                error_msg = f"{''.join(traceback.format_exception(exc_type, exc_value, exc_traceback))}"
                self.log_file.write(error_msg)
                self.log_file.flush()
                self.cleanup()
        except:  # noqa: E722, S110
            pass
        finally:
            self.original_excepthook(exc_type, exc_value, exc_traceback)


# see test_tracked.py for tests
def serialize_params_to_json(params: dict) -> dict:
    serialized_params = {}
    for key, value in params.items():
        # None and empty list are missing/empty values, skip them consistent with elsewhere in the code
        if value is None or (isinstance(value, list) and len(value) == 0):
            continue
        dtype, converted_value, _ = infer_convert_dtype_key_value(key, value, mute=True)
        # converted_value is not JSON if dtype is a SQLRecord or a list of SQLRecords
        # because we just the above function for features where we'd like to keep SQLRecords as they are
        # so, need to handle this here
        if (
            dtype == "?" or dtype.startswith("cat") or dtype.startswith("list[cat")
        ) and dtype not in {"cat ? str", "list[cat ? str]"}:
            if isinstance(value, SQLRecord):
                serialized_params[key] = (
                    f"{value.__class__.__get_name_with_module__()}[{value.uid}]"
                )
            elif dtype.startswith("list[cat"):
                items = list(value)
                if items and all(isinstance(item, SQLRecord) for item in items):
                    serialized_params[key] = [  # type: ignore
                        f"{item.__class__.__get_name_with_module__()}[{item.uid}]"
                        for item in items
                    ]
        else:
            serialized_params[key] = converted_value
        if key not in serialized_params:
            logger.warning(
                f"skipping param {key} with value {value} and dtype {dtype} not JSON serializable"
            )
            continue
        if is_sensitive_param_key(key) or is_sensitive_param_value(
            serialized_params[key]
        ):
            serialized_params[key] = REDACTED_SECRET_VALUE
    return serialized_params


class Context:
    """Run context.

    Is the book keeper for :func:`~lamindb.track` and :func:`~lamindb.finish`.
    """

    def __init__(self, uid: str | None = None, path: Path | None = None):
        self._uid: str | None = uid
        self._path: Path | None = path
        self._description: str | None = None
        self._version: str | None = None
        self._transform: Transform | None = None
        self._run: Run | None = None
        self._project: Project | None = None
        self._space: Space | None = None
        self._branch: Branch | None = None
        self._logging_message_track: str = ""
        self._logging_message_imports: str = ""
        self._stream_tracker: LogStreamTracker = LogStreamTracker()
        self._is_finish_retry: bool = False
        self._notebook_runner: str | None = None
        self._is_step_decorator_run: bool = False

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
    def space(self) -> Space | None:
        """The space in which artifacts, collections, transforms, and runs are saved during the run."""
        return self._space

    @property
    def branch(self) -> Branch | None:
        """The branch on which entities are created during the run."""
        return self._branch

    @property
    def run(self) -> Run | None:
        """Managed run of context."""
        return self._run

    def _track(
        self,
        transform: str | Transform | None = None,
        *,
        project: str | Project | None = None,
        space: str | Space | None = None,
        branch: str | Branch | None = None,
        plan: str | Artifact | None = None,
        features: dict | None = None,
        params: dict | None = None,
        new_run: bool | None = None,
        pypackages: bool | None = None,
        key: str | None = None,
        path: str | Path | None = None,
        source_code: str | None = None,
        kind: TransformKind | None = None,
        entrypoint: str | None = None,
        initiated_by_run: Run | str | None = None,
        stream_tracking: bool | None = None,
    ) -> None:
        """Track a run of a notebook or script.

        Populates the global run :class:`~lamindb.context` with :class:`~lamindb.Transform` & :class:`~lamindb.Run` objects and tracks the compute environment.

        Args:
            transform: A transform (stem) `uid` or object. If `None`, auto-creates a `transform` with its `uid`.
            project: A project or its `name` or `uid` for labeling entities created during the run.
            space: A restricted space or its `name` or `uid` in which to store entities created during the run.
                Default: the `"all"` space. Note that bionty entities ignore this setting and always get written to the `"all"` space.
                If you want to manually move entities to a different space, set the `.space` field (:doc:`docs:permissions`).
            branch: A branch (or its `name` or `uid`) on which to store records.
            plan: A plan, typically an agent plan. Pass an artifact (or its `key` or `uid`).
            features: A dictionary of features & values to track for the run.
            params: A dictionary of params & values to track for the run.
            new_run: If `False`, loads the latest run of transform
                (default notebook), if `True`, creates new run (default non-notebook).
            pypackages: If `True` or `None`, infers Python packages used in a notebook.
            key: Transform key.
            path: Filepath of a notebook or script.
            source_code: Source code.
            kind: Transform kind.
            entrypoint: Optional entrypoint name (e.g. function qualname) for the run.
            initiated_by_run: Optional parent run (or its `uid`) that triggered this run.
                If `None`, falls back to the `LAMIN_INITIATED_BY_RUN_UID` environment variable when set.
            stream_tracking: If set, override whether to capture stdout/stderr to run logs.
                Used by the flow/step decorator: flows get logs (`True`), steps do not (`False`).

        Examples
        --------

        To track the run of a notebook or script:

        .. literalinclude:: scripts/run_track_and_finish.py
            :language: python

        To ensure one version history across file renames::

            ln.track("Onv04I53OgtT")

        To track a project or an agent plan: pass a project/artifact to `ln.track()`, for example::

            ln.track(project="My project", plan="./plans/curate-dataset-x.md")

        Note that you have to create a project or save the agent plan in case it they don't yet exist::

            # create a project in Python
            ln.Project(name="My project").save()

            # create a project with the CLI
            lamin create project "My project"

            # save an agent plan with the CLI
            lamin save /path/to/.cursor/plans/curate-dataset-x.plan.md
            lamin save /path/to/.claude/plans/curate-dataset-x.md

        To sync code with a git repo, see: :ref:`sync-code-with-git`.

        To track parameters and features, see: :ref:`track-run-parameters`.

        To browse more examples, see: :doc:`/track`.
        """
        from lamindb.models import Artifact, Branch, Project, Space

        from .._finish import (
            save_context_core,
        )

        # similar logic here: https://github.com/laminlabs/lamindb/pull/2527
        if ln_setup.settings.instance.is_read_only_connection:
            logger.warning("skipping track(), connected in read-only mode")
            return None
        if project is None:
            project = os.environ.get("LAMIN_CURRENT_PROJECT")
        if project is not None:
            if isinstance(project, Project):
                assert project._state.adding is False, (  # noqa: S101
                    "Project must be saved before passing it to track()"
                )
                project_record = project
            else:
                project_record = Project.filter(
                    Q(name=project) | Q(uid=project)
                ).one_or_none()
                if project_record is None:
                    raise InvalidArgument(
                        f"Project '{project}' not found, either create it with `ln.Project(name='...').save()` or fix typos."
                    )
            self._project = project_record
        if space is not None:
            if isinstance(space, Space):
                assert space._state.adding is False, (  # noqa: S101
                    "Space must be saved before passing it to track()"
                )
                space_record = space
            else:
                space_record = Space.filter(Q(name=space) | Q(uid=space)).one_or_none()
                if space_record is None:
                    raise InvalidArgument(
                        f"Space '{space}', please check on the hub UI whether you have the correct `uid` or `name`."
                    )
            self._space = space_record
        if branch is not None:
            if isinstance(branch, Branch):
                assert branch._state.adding is False, (  # noqa: S101
                    "Branch must be saved before passing it to track()"
                )
                branch_record = branch
            else:
                branch_record = Branch.filter(
                    Q(name=branch) | Q(uid=branch)
                ).one_or_none()
                if branch_record is None:
                    raise InvalidArgument(
                        f"Space '{branch}', please check on the hub UI whether you have the correct `uid` or `name`."
                    )
            self._branch = branch_record
        plan_record: Artifact | None = None
        if plan is not None:
            if isinstance(plan, Artifact):
                assert plan._state.adding is False, (  # noqa: S101
                    "Plan artifact must be saved before passing it to track()"
                )
                plan_record = plan
            else:
                plan_record = Artifact.filter(Q(key=plan) | Q(uid=plan)).one_or_none()
                if plan_record is None:
                    raise InvalidArgument(
                        f"Plan artifact '{plan}' not found, either create it or use a valid key/uid."
                    )
        if initiated_by_run is None:
            initiated_by_run = os.environ.get("LAMIN_INITIATED_BY_RUN_UID")
        initiated_by_run_record: Run | None = None
        if initiated_by_run is not None:
            if isinstance(initiated_by_run, Run):
                assert initiated_by_run._state.adding is False, (  # noqa: S101
                    "initiated_by_run must be saved before passing it to track()"
                )
                initiated_by_run_record = initiated_by_run
            else:
                initiated_by_run_record = Run.filter(uid=initiated_by_run).one_or_none()
                if initiated_by_run_record is None:
                    raise InvalidArgument(
                        f"Run '{initiated_by_run}' not found, please pass a valid run uid."
                    )
        self._logging_message_track = ""
        self._logging_message_imports = ""
        self._is_step_decorator_run = (
            entrypoint is not None and stream_tracking is False
        )
        if transform is not None and isinstance(transform, str):
            self.uid = transform
            transform = None
            uid_was_none = False
        else:
            uid_was_none = True
        self._path = None
        cli_call = get_cli_call()
        if transform is None:
            description = None
            if source_code is not None:
                transform_kind = kind if kind is not None else "function"
                assert key is not None, (
                    "`key` cannot be `None` when `source_code` is passed to `track()`."
                )
                assert path is None, (
                    "`path` cannot be passed when `source_code` is passed to `track()`."
                )
            else:
                if is_run_from_ipython or _is_running_in_marimo():
                    self._path, description = self._track_notebook(
                        path_str=path, pypackages=pypackages
                    )
                    transform_kind = "notebook"
                else:
                    (
                        self._path,
                        transform_kind,
                        _,
                        _,
                        key_from_module,
                    ) = detect_and_process_source_code_file(
                        path=path, infer_reference=False
                    )
                    if key is None and key_from_module is not None:
                        key = key_from_module
            if description is None:
                description = self._description
            if description is None and cli_call is not None:
                description = f"CLI: {cli_call[0]}"
            self._create_or_load_transform(
                description=description,
                transform_kind=transform_kind,
                key=key,
                source_code=source_code,
            )
        else:
            if transform.kind in {"notebook", "script"}:
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
                self._logging_message_track += (
                    f"created Transform('{transform.uid}', key='{transform.key}')"
                )
                transform_exists = transform
            else:
                self._logging_message_track += (
                    f"loaded Transform('{transform.uid}', key='{transform.key}')"
                )
            self._transform = transform_exists

        if new_run is None:  # for notebooks, default to loading latest runs
            new_run = (
                False
                if (
                    self._transform.kind == "notebook"
                    and self._notebook_runner != "nbconvert"
                )
                else True
            )  # type: ignore

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
                run._status_code = -2  # re-started
                if plan_record is not None:
                    run.plan = plan_record
                    run.save()
                entrypoint_str = (
                    f", entrypoint='{entrypoint}'" if entrypoint is not None else ""
                )
                self._logging_message_track += f", re-started Run('{run.uid}'{entrypoint_str}) at {format_field_value(run.started_at)}"

        if run is None:  # create new run with auto-populated started_at timestamp
            run = Run(
                transform=self._transform,
                status="started",
                initiated_by_run=initiated_by_run_record,
                entrypoint=entrypoint,
                plan=plan_record,
            )
            entrypoint_str = (
                f", entrypoint='{entrypoint}'" if entrypoint is not None else ""
            )
            # run.started_at is set on insert by the database (db_default=Now()), independently of the Python code
            # hence we log datetime.now(timezone.utc) instead of run.started_at
            self._logging_message_track += f", started new Run('{run.uid}'{entrypoint_str}) at {format_field_value(datetime.now(timezone.utc))}"
        # can only determine at ln.finish() if run was consecutive in
        # interactive session, otherwise, is consecutive
        run.is_consecutive = True if is_run_from_ipython else None
        if params is not None:
            run.params = serialize_params_to_json(params)
            self._logging_message_track += "\n→ params: " + ", ".join(
                f"{key}={value!r}" for key, value in run.params.items()
            )
        if cli_call is not None:
            _, cli_args = cli_call
            logger.important(f"script invoked with: {cli_args}")
            run.cli_args = cli_args
        run.save()  # need to save now
        if features is not None:
            run.features.add_values(features)
            self._logging_message_track += "\n→ features: " + ", ".join(
                f"{key}={value!r}" for key, value in features.items()
            )
        self._run = run
        track_python_environment(run)
        if self.project is not None:
            # to update a potential project link
            # is only necessary if transform is loaded rather than newly created
            # can be optimized by checking whether the transform is loaded, but it typically is
            self.transform.save()
        log_to_file = None
        if log_to_file is None:
            if stream_tracking is not None:
                log_to_file = stream_tracking
            else:
                # Script runs get stream tracking; decorator-based runs only when
                # stream_tracking is passed (flow=True from decorator).
                log_to_file = self.transform.kind == "script"
        if log_to_file:
            self._stream_tracker.start(run)
        logger.important(self._logging_message_track)
        if self._logging_message_imports:
            logger.important(self._logging_message_imports)
        if uid_was_none and self._path is not None:
            # Flow/step decorators set run.entrypoint. Show this recommendation only
            # for flows (`stream_tracking=True`) and suppress it for steps.
            if entrypoint is not None:
                if stream_tracking:
                    logger.important_hint(
                        f'recommendation: to identify the script across renames, pass the uid: @ln.flow(uid="{self.transform.uid[:-4]}")'
                    )
            else:
                notebook_or_script = (
                    "notebook" if self._transform.kind == "notebook" else "script"
                )
                r_or_python = "." if self._path.suffix in {".py", ".ipynb"} else "$"
                project_str = (
                    f', project="{project if isinstance(project, str) else project.name}"'
                    if project is not None
                    else ""
                )
                space_str = (
                    f', space="{space if isinstance(space, str) else space.name}"'
                    if space is not None
                    else ""
                )
                plan_str = (
                    f', plan="{plan if isinstance(plan, str) else plan.key}"'
                    if plan is not None
                    else ""
                )
                params_str = (
                    ", params={...}" if params is not None else ""
                )  # do not put the values because typically parameterized by user
                kwargs_str = f"{project_str}{space_str}{plan_str}{params_str}"
                logger.important_hint(
                    f'recommendation: to identify the {notebook_or_script} across renames, pass the uid: ln{r_or_python}track("{self.transform.uid[:-4]}"{kwargs_str})'
                )
        if (
            self.transform.kind == "script"
            and self._path is not None
            and not self._is_step_decorator_run
        ):
            save_context_core(
                run=run,
                transform=self.transform,
                filepath=self._path,
                message_prefix="monitor at",
            )

    def _track_notebook(
        self,
        *,
        path_str: str | Path | None,
        pypackages: bool | None = None,
    ) -> tuple[Path, str | None]:
        if path_str is None:
            path, self._notebook_runner = get_notebook_path()
        else:
            path = Path(path_str)

        if pypackages is None:
            pypackages = True
        description = None

        if self._notebook_runner == "marimo":
            import re

            source = path.read_text(encoding="utf-8")
            auto_download_ipynb_re = re.compile(
                r"""
                            auto_download
                            \s*=\s*
                            \[
                            [^\]]*
                            ["']ipynb["']
                            [^\]]*
                            \]
                            """,
                re.VERBOSE,
            )
            if not auto_download_ipynb_re.search(source):
                raise SystemExit(
                    "Tracking marimo run reports requires auto-export of ipynb files.\n"
                    "In the top left corner of the UI, go to Settings → Exporting outputs and then select `ipynb`."
                )
            return path, description
        if path.suffix == ".ipynb" and path.stem.startswith("Untitled"):
            raise RuntimeError(
                "Your notebook file name is 'Untitled.ipynb', please rename it before tracking. You might have to re-start your notebook kernel."
            )
        path_str = path.as_posix()
        if path_str.startswith("/fileId="):
            logger.warning("tracking on Google Colab is experimental")
            path_str = get_notebook_key_colab()
            path = Path(path_str)
        else:
            from nbproject.dev import read_notebook
            from nbproject.dev._meta_live import get_title
            from nbproject.dev._pypackage import infer_pypackages

            try:
                nb = read_notebook(path_str)

                nbproject_title = get_title(nb)
                if nbproject_title is not None:
                    description = nbproject_title

                if pypackages:
                    self._logging_message_imports += (
                        "notebook imports:"
                        f" {pretty_pypackages(infer_pypackages(nb, pin_versions=True))}"
                    )
            except Exception:
                logger.debug("reading the notebook file failed")
                pass
        return path, description

    def _create_or_load_transform(
        self,
        *,
        description: str | None = None,
        transform_ref: str | None = None,
        transform_ref_type: str | None = None,
        transform_kind: TransformKind = None,
        key: str | None = None,
        source_code: str | None = None,
    ):
        transform, logging_message = Transform._create_or_load_from_source(
            path=self._path,
            description=description,
            transform_ref=transform_ref,
            transform_ref_type=transform_ref_type,
            transform_kind=transform_kind,
            key=key,
            source_code=source_code,
            uid=self.uid,
            version=self.version,
            notebook_runner=self._notebook_runner,
        )
        self._transform = transform
        self._uid = transform.uid
        self._logging_message_track += logging_message

    def _finish(self, ignore_non_consecutive: None | bool = None) -> None:
        """Finish the run of a notebook or script.

        - writes a timestamp: `run.finished_at`
        - saves the source code if it is not yet saved: `transform.source_code`
        - saves a run report: `run.report`

        When called in a notebook, will prompt to save the notebook in your editor.

        In a Jupyter notebook, call `ln.finish()` in its own cell as the output of the cell in which
        `ln.finish()` is called is stripped from the run report.

        Args:
            ignore_non_consecutive: Whether to ignore if a notebook was non-consecutively executed.

        See Also:
            `lamin save script.py` or `lamin save notebook.ipynb` → `docs </cli#lamin-save>`__

        Examples
        --------

        See :doc:`/track`.

        """
        from .._finish import save_context_core, save_run_logs

        if self.run is None:
            raise TrackNotCalled("Please run `ln.track()` before `ln.finish()`")
        if self._path is None:
            if self.run.transform.kind in {"script", "notebook"}:
                raise ValueError(
                    "Transform type is not allowed to be 'script' or 'notebook' because `context._path` is `None`."
                )
            self.run.finished_at = datetime.now(timezone.utc)
            self.run.save()
            # reset context so the next _track() starts clean (e.g. from decorator)
            self._uid = None
            self._run = None
            self._transform = None
            self._version = None
            self._description = None
            self._is_step_decorator_run = False
            return None
        self.run._status_code = 0
        if self.transform.kind == "notebook":
            return_code = save_context_core(
                run=self.run,
                transform=self.run.transform,
                filepath=self._path,
                finished_at=True,
                ignore_non_consecutive=ignore_non_consecutive,
                is_retry=self._is_finish_retry,
                notebook_runner=self._notebook_runner,
            )
            if return_code == "retry":
                self._is_finish_retry = True
                return None
        else:
            self.run.finished_at = datetime.now(timezone.utc)
            self.run.save()  # persist finished_at (save_run_logs only saves when log file exists)
            if ln_setup.settings.instance.is_on_hub and not self._is_step_decorator_run:
                instance_slug = ln_setup.settings.instance.slug
                ui_url = ln_setup.settings.instance.ui_url
                logger.important(f"go to: {ui_url}/{instance_slug}/run/{self.run.uid}")
            save_run_logs(self.run, save_run=True)
            self._stream_tracker.finish()
        # reset the context attributes so that somebody who runs `track()` after finish
        # starts fresh
        self._uid = None
        self._run = None
        self._transform = None
        self._version = None
        self._description = None
        self._is_step_decorator_run = False


context: Context = Context()
