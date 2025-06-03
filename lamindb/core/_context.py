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
    InvalidArgument,
    TrackNotCalled,
    UpdateContext,
)
from ..models._is_versioned import bump_version as bump_version_function
from ..models._is_versioned import (
    increment_base62,
)
from ._sync_git import get_transform_reference_from_git_repo
from ._track_environment import track_environment

if TYPE_CHECKING:
    from lamindb_setup.core.types import UPathStr

    from lamindb.base.types import TransformType
    from lamindb.models import Branch, Project, Space

is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)

msg_path_failed = "failed to infer notebook path.\nfix: pass `path` to `ln.track()`"


def get_uid_ext(version: str) -> str:
    from lamin_utils._base62 import encodebytes

    # merely zero-padding the nbproject version such that the base62 encoding is
    # at least 4 characters long doesn't yields sufficiently diverse hashes and
    # leads to collisions; it'd be nice because the uid_ext would be ordered
    return encodebytes(hashlib.md5(version.encode()).digest())[:4]  # noqa: S324


def get_notebook_path() -> tuple[Path, str]:
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

    Is the book keeper for :meth:`~lamindb.core.Context.track`.
    """

    def __init__(self):
        self._uid: str | None = None
        self._description: str | None = None
        self._version: str | None = None
        self._transform: Transform | None = None
        self._run: Run | None = None
        self._path: Path | None = None
        """A local path to the script or notebook that's running."""
        self._project: Project | None = None
        self._space: Space | None = None
        self._branch: Branch | None = None
        self._logging_message_track: str = ""
        self._logging_message_imports: str = ""
        self._stream_tracker: LogStreamTracker = LogStreamTracker()
        self._is_finish_retry: bool = False
        self._notebook_runner: str | None = None

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
        params: dict | None = None,
        new_run: bool | None = None,
        path: str | None = None,
    ) -> None:
        """Track a run of your notebook or script.

        Populates the global run :class:`~lamindb.context` by managing `Transform` & `Run` records and caching the compute environment.

        If :attr:`~lamindb.core.Settings.sync_git_repo` is set, checks whether a script-like transform exists in a git repository and links it.

        Args:
            transform: A transform (stem) `uid` (or record). If `None`, auto-creates a `transform` with its `uid`.
            project: A project (or its `name` or `uid`) for labeling entities.
            space: A restricted space (or its `name` or `uid`) in which to store artifacts, collections, transforms, and runs.
                Default: the `"All"` space.
                If you want to manually move entities to a different space, set the `.space` field (:doc:`docs:access`).
            branch: A branch (or its `name` or `uid`) on which to store records.
            params: A dictionary of parameters to track for the run.
            new_run: If `False`, loads the latest run of transform
                (default notebook), if `True`, creates new run (default non-notebook).
            path: Filepath of notebook or script. Only needed if it can't be
                automatically detected.

        Examples:

            To track the run of a notebook or script, call::

                ln.track()
                #> → created Transform('Onv04I53OgtT0000'), started new Run('dpSfd7Ds...') at 2025-04-25 11:00:03 UTC
                #> • recommendation: to identify the notebook across renames, pass the uid: ln.track("Onv04I53OgtT")

            Ensure one version history across file renames::

                ln.track("Onv04I53OgtT")
                #> → created Transform('Onv04I53OgtT0000'), started new Run('dpSfd7Ds...') at 2025-04-25 11:00:03 UTC

            More examples: :doc:`/track`
        """
        from lamindb.models import Branch, Project, Space

        instance_settings = ln_setup.settings.instance
        # similar logic here: https://github.com/laminlabs/lamindb/pull/2527
        # TODO: refactor upon new access management
        if instance_settings.dialect == "postgresql" and "read" in instance_settings.db:
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
        self._logging_message_track = ""
        self._logging_message_imports = ""
        if transform is not None and isinstance(transform, str):
            self.uid = transform
            transform = None
            uid_was_none = False
        else:
            uid_was_none = True
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
            new_run = (
                False
                if (
                    self._transform.type == "notebook"
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
            run.features.add_values(params)
            self._logging_message_track += "\n→ params: " + ", ".join(
                f"{key}={value}" for key, value in params.items()
            )
        self._run = run
        track_environment(run)
        if self.project is not None:
            # to update a potential project link
            # is only necessary if transform is loaded rather than newly created
            # can be optimized by checking whether the transform is loaded, but it typically is
            self.transform.save()
        log_to_file = None
        if log_to_file is None:
            log_to_file = self.transform.type != "notebook"
        if log_to_file:
            self._stream_tracker.start(run)
        logger.important(self._logging_message_track)
        if self._logging_message_imports:
            logger.important(self._logging_message_imports)
        if uid_was_none:
            notebook_or_script = (
                "notebook" if self._transform.type == "notebook" else "script"
            )
            r_or_python = "."
            if self._path is not None:
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
            params_str = (
                ", params={...}" if params is not None else ""
            )  # do not put the values because typically parameterized by user
            kwargs_str = f"{project_str}{space_str}{params_str}"
            logger.important_hint(
                f'recommendation: to identify the {notebook_or_script} across renames, pass the uid: ln{r_or_python}track("{self.transform.uid[:-4]}"{kwargs_str})'
            )

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
            path, self._notebook_runner = get_notebook_path()
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
            from nbproject.dev import read_notebook
            from nbproject.dev._meta_live import get_title
            from nbproject.dev._pypackage import infer_pypackages

            try:
                nb = read_notebook(path_str)

                nbproject_title = get_title(nb)
                if nbproject_title is not None:
                    description = nbproject_title

                self._logging_message_imports += (
                    "notebook imports:"
                    f" {pretty_pypackages(infer_pypackages(nb, pin_versions=True))}"
                )
            except Exception:
                logger.debug("reading the notebook file failed")
                pass
        return path, description

    def _process_aux_transform(
        self,
        aux_transform: Transform,
        transform_hash: str,
    ) -> tuple[str, Transform | None, str]:
        # first part of the if condition: no version bump, second part: version bump
        message = ""
        if (
            # if a user hasn't yet saved the transform source code AND is the same user
            (
                aux_transform.source_code is None
                and aux_transform.created_by_id == ln_setup.settings.user.id
            )
            # if the transform source code is unchanged
            # if aux_transform.type == "notebook", we anticipate the user makes changes to the notebook source code
            # in an interactive session, hence we *pro-actively bump* the version number by setting `revises` / 'nbconvert' execution is NOT interactive
            # in the second part of the if condition even though the source code is unchanged at point of running track()
            or (
                aux_transform.hash == transform_hash
                and (
                    aux_transform.type != "notebook"
                    or self._notebook_runner == "nbconvert"
                )
            )
        ):
            uid = aux_transform.uid
            return uid, aux_transform, message
        else:
            uid = f"{aux_transform.uid[:-4]}{increment_base62(aux_transform.uid[-4:])}"
            message = (
                f"found {aux_transform.type} {aux_transform.key}, making new version"
            )
            if (
                aux_transform.hash == transform_hash
                and aux_transform.type == "notebook"
            ):
                message += " -- anticipating changes"
            elif aux_transform.hash != transform_hash:
                message += (
                    ""  # could log "source code changed", but this seems too much
                )
            elif aux_transform.created_by_id != ln_setup.settings.user.id:
                message += (
                    f" -- {aux_transform.created_by.handle} already works on this draft"
                )
            return uid, None, message

    def _create_or_load_transform(
        self,
        *,
        description: str,
        transform_ref: str | None = None,
        transform_ref_type: str | None = None,
        transform_type: TransformType = None,
    ):
        from .._finish import notebook_to_script

        if not self._path.suffix == ".ipynb":
            transform_hash, _ = hash_file(self._path)
        else:
            # need to convert to stripped py:percent format for hashing
            source_code_path = ln_setup.settings.cache_dir / self._path.name.replace(
                ".ipynb", ".py"
            )
            notebook_to_script(description, self._path, source_code_path)
            transform_hash, _ = hash_file(source_code_path)
        # see whether we find a transform with the exact same hash
        aux_transform = Transform.filter(hash=transform_hash).one_or_none()
        # if the user did not pass a uid and there is no matching aux_transform
        # need to search for the transform based on the filename
        if self.uid is None and aux_transform is None:

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
            if len(transforms) != 0:
                message = ""
                found_key = False
                for aux_transform in transforms:
                    if aux_transform.key in self._path.as_posix():
                        key = aux_transform.key
                        uid, target_transform, message = self._process_aux_transform(
                            aux_transform, transform_hash
                        )
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
        elif self.uid is not None and len(self.uid) == 16:
            transform = Transform.filter(uid=self.uid).one_or_none()
            if transform is not None:
                if transform.key not in self._path.as_posix():
                    n_parts = len(Path(transform.key).parts)
                    (
                        Path(*self._path.parts[-n_parts:]).as_posix()
                        if n_parts > 0
                        else ""
                    )
                    key = self._path.name
                else:
                    key = transform.key  # type: ignore
            else:
                key = self._path.name
        else:
            if self.uid is not None:
                # the case with length 16 is covered above
                if not len(self.uid) == 12:
                    raise InvalidArgument(
                        f'Please pass an auto-generated uid instead of "{self.uid}". Resolve by running: ln.track("{base62_12()}")'
                    )
                aux_transform = (
                    Transform.filter(uid__startswith=self.uid)
                    .order_by("-created_at")
                    .first()
                )
            else:
                # deal with a hash-based match
                # the user might have a made a copy of the notebook or script
                # and actually wants to create a new transform
                if aux_transform is not None and not aux_transform.key.endswith(
                    self._path.name
                ):
                    prompt = f"Found transform with same hash but different key: {aux_transform.key}. Did you rename your {transform_type} to {self._path.name} (1) or intentionally made a copy (2)?"
                    response = (
                        "1" if os.getenv("LAMIN_TESTING") == "true" else input(prompt)
                    )
                    assert response in {"1", "2"}, (  # noqa: S101
                        f"Please respond with either 1 or 2, not {response}"
                    )
                    if response == "2":
                        transform_hash = None  # make a new transform
            if aux_transform is not None:
                if aux_transform.key.endswith(self._path.name):
                    key = aux_transform.key
                else:
                    key = "/".join(
                        aux_transform.key.split("/")[:-1] + [self._path.name]
                    )
                uid, target_transform, message = self._process_aux_transform(
                    aux_transform, transform_hash
                )
                if message != "":
                    logger.important(message)
            else:
                uid = f"{self.uid}0000" if self.uid is not None else None
                target_transform = None
                key = self._path.name
            self.uid, transform = uid, target_transform
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
            if self.uid is not None and len(self.uid) == 16:
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
            transform = Transform(  # type: ignore
                uid=self.uid,
                version=self.version,
                description=description,
                key=key,
                reference=transform_ref,
                reference_type=transform_ref_type,
                type=transform_type,
            ).save()
            self._logging_message_track += f"created Transform('{transform.uid}')"
        else:
            uid = transform.uid
            # transform was already saved via `finish()`
            transform_was_saved = transform.source_code is not None
            # check whether the transform.key is consistent
            if transform.key != key:
                self._logging_message_track += (
                    f"renaming transform {transform.key} to {key}"
                )
                transform.key = key
                transform.save()
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
                    f'{transform.created_by.name} ({transform.created_by.handle}) already works on this draft {transform.type}.\n\nPlease create a revision via `ln.track("{uid[:-4]}{increment_base62(uid[-4:])}")` or a new transform with a *different* key and `ln.track("{ids.base62_12()}0000")`.'
                )
            # check whether transform source code was already saved
            if transform_was_saved:
                bump_revision = False
                if (
                    transform.type == "notebook"
                    and self._notebook_runner != "nbconvert"
                ):
                    # we anticipate the user makes changes to the notebook source code
                    # in an interactive session, hence we pro-actively bump the version number
                    bump_revision = True
                else:
                    if transform_hash != transform.hash:
                        bump_revision = True
                    else:
                        self._logging_message_track += (
                            f"loaded Transform('{transform.uid}')"
                        )
                if bump_revision:
                    change_type = (
                        "re-running notebook with already-saved source code"
                        if (
                            transform.type == "notebook"
                            and self._notebook_runner != "nbconvert"
                        )
                        else "source code changed"
                    )
                    raise UpdateContext(
                        f'✗ {change_type}, please update the `uid` argument in `track()` to "{uid[:-4]}{increment_base62(uid[-4:])}"'
                    )
            else:
                self._logging_message_track += f"loaded Transform('{transform.uid}')"
        self._transform = transform

    def _finish(self, ignore_non_consecutive: None | bool = None) -> None:
        """Finish the run and write a run report.

        - writes a timestamp: `run.finished_at`
        - saves the source code if it is not yet saved: `transform.source_code`
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
            notebook_runner=self._notebook_runner,
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

    @deprecated("ln.track()")
    def track(self, *args, **kwargs):
        return self._track(*args, **kwargs)

    @deprecated("ln.finish()")
    def finish(self, *args, **kwargs):
        return self._finish(*args, **kwargs)


context: Context = Context()
