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
from lnschema_core.models import Param, ParamValue, RunParamValue
from lnschema_core.types import TransformType
from lnschema_core.users import current_user_id

from lamindb.core._transform_settings import transform as transform_settings

from ._settings import settings
from ._sync_git import get_transform_reference_from_git_repo
from .exceptions import (
    MissingTransformSettings,
    NotebookNotSavedError,
    NoTitleError,
    UpdateTransformSettings,
)
from .versioning import bump_version as bump_version_function

if TYPE_CHECKING:
    from lamindb_setup.core.types import UPathStr

is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)

msg_path_failed = (
    "failed to infer notebook path.\nfix: either track manually via"
    " `ln.track(transform=ln.Transform(name='My notebook'))` or pass"
    " `path` to ln.track()"
)


def get_uid_ext(version: str) -> str:
    from lamin_utils._base62 import encodebytes

    # merely zero-padding the nbproject version such that the base62 encoding is
    # at least 4 characters long doesn't yields sufficiently diverse hashes and
    # leads to collisions; it'd be nice because the uid_ext would be ordered
    return encodebytes(hashlib.md5(version.encode()).digest())[:4]


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
            new_version = bump_version_function(version, behavior="prompt")
        else:
            new_version = response
        updated = new_version != version
    if updated:
        new_metadata = (
            f'ln.settings.transform.stem_uid = "{new_stem_uid}"\nln.settings.transform.version ='
            f' "{new_version}"\n'
        )
        raise UpdateTransformSettings(
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


MESSAGE = """To track this {transform_type}, set

ln.settings.transform.stem_uid = "{stem_uid}"
ln.settings.transform.version = "{version}"
"""

MESSAGE_UPDATE = """You updated your {transform_type}.

If this is a minor update, bump your version from {old_version} to:

ln.settings.transform.version = "{new_version_minor_bump}"

If this is a major update, bump it to:

ln.settings.transform.version = "{new_version_major_bump}"

If this is a new {transform_type}, set:

ln.settings.transform.stem_uid = "{new_stem_uid}"
ln.settings.transform.version = "1"

"""


def raise_transform_settings_error_needs_update(old_version: str) -> None:
    from lnschema_core.ids import base62_12

    transform_type = "notebook" if is_run_from_ipython else "script"
    new_stem_uid = base62_12()

    raise UpdateTransformSettings(
        MESSAGE_UPDATE.format(
            transform_type=transform_type,
            new_stem_uid=new_stem_uid,
            old_version=old_version,
            new_version_major_bump=bump_version_function(
                old_version, bump_type="major", behavior="ignore"
            ),
            new_version_minor_bump=bump_version_function(
                old_version, bump_type="minor", behavior="ignore"
            ),
        )
    )


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
    raise MissingTransformSettings(
        MESSAGE.format(
            transform_type=transform_type, stem_uid=stem_uid, version=version
        )
    )


def pretty_pypackages(dependencies: dict) -> str:
    deps_list = []
    for pkg, ver in dependencies.items():
        if ver != "":
            deps_list.append(pkg + f"=={ver}")
        else:
            deps_list.append(pkg)
    deps_list.sort()
    return " ".join(deps_list)


def parse_and_link_params(run: Run, params: dict) -> None:
    param_values = []
    for key, value in params.items():
        param = Param.filter(name=key).one_or_none()
        if param is None:
            dtype = type(value).__name__
            logger.warning(
                f"param '{key}' does not yet exist, creating it with dtype '{dtype}'"
            )
            param = Param(name=key, dtype=dtype).save()
        param_value, _ = ParamValue.objects.get_or_create(param=param, value=value)
        param_values.append(param_value)
    if param_values:
        links = [
            RunParamValue(run_id=run.id, paramvalue_id=param_value.id)
            for param_value in param_values
        ]
        RunParamValue.objects.bulk_create(links)


class run_context:
    """Global run context."""

    transform: Transform | None = None
    """Current transform."""
    run: Run | None = None
    """Current run."""
    path: Path | None = None
    """A local path to the script that's running."""

    @classmethod
    def _track(
        cls,
        *,
        params: dict | None = None,
        transform: Transform | None = None,
        new_run: bool | None = None,
        path: str | None = None,
    ) -> Run:
        """Track notebook or script run.

        Creates or loads a global :class:`~lamindb.Run` that enables data
        lineage tracking. You can find it in :class:`~lamindb.core.run_context`.

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
            >>> ln.track()

            If you'd like to track an abstract pipeline run, pass a
            :class:`~lamindb.Transform` object of ``type`` ``"pipeline"``:

            >>> ln.Transform(name="Cell Ranger", version="2", type="pipeline").save()
            >>> transform = ln.Transform.filter(name="Cell Ranger", version="2").one()
            >>> ln.track(transform=transform)
        """
        cls.path = None
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
                    key, name = cls._track_notebook(path=path)
                    transform_type = TransformType.notebook
                    transform_ref = None
                    transform_ref_type = None
                else:
                    (name, key, transform_ref, transform_ref_type) = cls._track_script(
                        path=path
                    )
                    transform_type = TransformType.script
                cls._create_or_load_transform(
                    stem_uid=stem_uid,
                    version=version,
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
                raise_transform_settings_error()
        else:
            if transform.type in {"notebook", "script"}:
                raise ValueError(
                    "Use ln.track() without passing transform in a notebook or script"
                    " - metadata is automatically parsed"
                )
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
                run.started_at = datetime.now(timezone.utc)  # update run time
                logger.important(f"loaded: {run}")

        if run is None:  # create new run
            run = Run(
                transform=cls.transform,
                params=params,
            )
            logger.important(f"saved: {run}")
        # can only determine at ln.finish() if run was consecutive in
        # interactive session, otherwise, is consecutive
        run.is_consecutive = True if is_run_from_ipython else None
        # need to save in all cases
        run.save()
        if params is not None:
            parse_and_link_params(run, params)
        cls.run = run

        from ._track_environment import track_environment

        track_environment(run)
        return run

    @classmethod
    def _track_script(
        cls,
        *,
        path: UPathStr | None,
    ) -> tuple[str, str, str, str]:
        if path is None:
            import inspect

            frame = inspect.stack()[2]
            module = inspect.getmodule(frame[0])
            cls.path = Path(module.__file__)
        else:
            cls.path = Path(path)
        name = cls.path.name
        key = name
        reference = None
        reference_type = None
        if settings.sync_git_repo is not None:
            reference = get_transform_reference_from_git_repo(cls.path)
            reference_type = "url"
        return name, key, reference, reference_type

    @classmethod
    def _track_notebook(
        cls,
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
                from nbproject.dev._pypackage import infer_pypackages

                nb = nbproject.dev.read_notebook(path_str)
                logger.important(
                    "notebook imports:"
                    f" {pretty_pypackages(infer_pypackages(nb, pin_versions=True))}"
                )
            except Exception:
                logger.debug("inferring imported packages failed")
                pass
        cls.path = Path(path_str)
        return key, name

    @classmethod
    def _create_or_load_transform(
        cls,
        *,
        stem_uid: str,
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
            logger.important(f"saved: {transform}")
        else:
            # check whether there was an update to the transform, like
            # renaming the file or updating the title
            if transform.name != name or transform.key != key:
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
                    transform.key = key
                    transform.save()
                    logger.important(f"updated: {transform}")
            # check whether transform source code was already saved
            if transform.source_code_id is not None:
                response = None
                if is_run_from_ipython:
                    if os.getenv("LAMIN_TESTING") is None:
                        response = input(
                            "You already saved source code for this notebook."
                            " Auto-bump the version before a new run? (y/n)"
                        )
                    else:
                        response = "y"
                else:
                    hash, _ = hash_file(cls.path)  # ignore hash_type for now
                    if hash != transform.source_code.hash:
                        # only if hashes don't match, we need user input
                        if os.getenv("LAMIN_TESTING") is None:
                            response = input(
                                "You already saved source code for this script and meanwhile modified it without bumping a version."
                                " Auto-bump the version before a new run? (y/n)"
                            )
                        else:
                            response = "y"
                    else:
                        logger.important(f"loaded: {transform}")
                if response is not None:
                    # if a script is re-run and hashes match, we don't need user input
                    if response == "y":
                        update_stem_uid_or_version(stem_uid, version, bump_version=True)
                    else:
                        # the user didn't agree to auto-bump, hence treat manually
                        raise_transform_settings_error_needs_update(
                            old_version=transform.version
                        )
            else:
                logger.important(f"loaded: {transform}")
        cls.transform = transform
