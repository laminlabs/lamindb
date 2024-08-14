from __future__ import annotations

import os
import shutil
import subprocess
from datetime import datetime, timezone
from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
from lamin_utils import logger
from lamindb_setup.core.hashing import hash_file
from lnschema_core.types import TransformType

from .core._run_context import context, is_run_from_ipython

if TYPE_CHECKING:
    from pathlib import Path

    from lnschema_core import Run, Transform

    from ._query_set import QuerySet


class TrackNotCalled(SystemExit):
    pass


class NotebookNotSaved(SystemExit):
    pass


def get_seconds_since_modified(filepath) -> float:
    return datetime.now().timestamp() - filepath.stat().st_mtime


def finish() -> None:
    """Mark a tracked run as finished.

    Saves source code and, for notebooks, a run report to your default storage location.
    """
    if context.run is None:
        raise TrackNotCalled("Please run `ln.track()` before `ln.finish()`")
    if context._path is None:
        if context.run.transform.type in {"script", "notebook"}:
            raise ValueError(
                f"Transform type is not allowed to be 'script' or 'notebook' but is {context.run.transform.type}."
            )
        context.run.finished_at = datetime.now(timezone.utc)
        context.run.save()
        # nothing else to do
        return None
    if is_run_from_ipython:  # notebooks
        if (
            get_seconds_since_modified(context._path) > 3
            and os.getenv("LAMIN_TESTING") is None
        ):
            raise NotebookNotSaved(
                "Please save the notebook in your editor right before running `ln.finish()`"
            )
    save_run_context_core(
        run=context.run,
        transform=context.run.transform,
        filepath=context._path,
        finished_at=True,
    )


def save_run_context_core(
    *,
    run: Run,
    transform: Transform,
    filepath: Path,
    transform_family: QuerySet | None = None,
    finished_at: bool = False,
    from_cli: bool = False,
) -> str | None:
    import lamindb as ln

    ln.settings.verbosity = "success"

    # for scripts, things are easy
    is_consecutive = True
    is_notebook = transform.type == "notebook"
    _source_code_artifact_path = filepath
    # for notebooks, we need more work
    if is_notebook:
        try:
            import nbstripout
            from nbproject.dev import (
                check_consecutiveness,
                read_notebook,
            )
        except ImportError:
            logger.error(
                "install nbproject & nbstripout: pip install nbproject nbstripout"
            )
            return None
        notebook_content = read_notebook(filepath)  # type: ignore
        is_consecutive = check_consecutiveness(
            notebook_content, calling_statement="ln.finish()"
        )
        if not is_consecutive:
            msg = "   Do you still want to proceed with finishing? (y/n) "
            if os.getenv("LAMIN_TESTING") is None:
                response = input(msg)
            else:
                response = "n"
            if response != "y":
                return "aborted-non-consecutive"
        # convert the notebook file to html
        # log_level is set to 40 to silence the nbconvert logging
        subprocess.run(
            [
                "jupyter",
                "nbconvert",
                "--to",
                "html",
                filepath.as_posix(),
                "--Application.log_level=40",
            ],
            check=True,
        )
        # move the temporary file into the cache dir in case it's accidentally
        # in an existing storage location -> we want to move associated
        # artifacts into default storage and not register them in an existing
        # location
        report_path_orig = filepath.with_suffix(".html")  # current location
        report_path = ln_setup.settings.storage.cache_dir / report_path_orig.name
        # don't use Path.rename here because of cross-device link error
        # https://laminlabs.slack.com/archives/C04A0RMA0SC/p1710259102686969
        shutil.move(
            report_path_orig,  # type: ignore
            report_path,
        )
        # strip the output from the notebook to create the source code file
        # first, copy the notebook file to a temporary file in the cache
        _source_code_artifact_path = ln_setup.settings.storage.cache_dir / filepath.name
        shutil.copy2(filepath, _source_code_artifact_path)  # copy
        subprocess.run(
            [
                "nbstripout",
                _source_code_artifact_path,
                "--extra-keys",
                "metadata.version metadata.kernelspec metadata.language_info metadata.pygments_lexer metadata.name metadata.file_extension",
            ],
            check=True,
        )
    # find initial versions of source codes and html reports
    prev_report = None
    prev_source = None
    if transform_family is None:
        transform_family = transform.versions
    if len(transform_family) > 0:
        for prev_transform in transform_family.order_by("-created_at"):
            if (
                prev_transform.latest_run is not None
                and prev_transform.latest_run.report_id is not None
            ):
                prev_report = prev_transform.latest_run.report
            if prev_transform._source_code_artifact_id is not None:
                prev_source = prev_transform._source_code_artifact
    ln.settings.creation.artifact_silence_missing_run_warning = True

    # track source code
    if transform._source_code_artifact_id is not None:
        # check if the hash of the transform source code matches
        # (for scripts, we already run the same logic in track() - we can deduplicate the call at some point)
        hash, _ = hash_file(_source_code_artifact_path)  # ignore hash_type for now
        if hash != transform._source_code_artifact.hash:
            if os.getenv("LAMIN_TESTING") is None:
                # in test, auto-confirm overwrite
                response = input(
                    f"You are about to replace (overwrite) existing source code (hash '{transform._source_code_artifact.hash}') for transform version"
                    f" '{transform.version}'. Proceed? (y/n)"
                )
            else:
                response = "y"
            if response == "y":
                transform._source_code_artifact.replace(_source_code_artifact_path)
                transform._source_code_artifact.save(upload=True)
                logger.success(
                    f"replaced transform._source_code_artifact: {transform._source_code_artifact}"
                )
            else:
                logger.warning("Please re-run `ln.track()` to make a new version")
                return "rerun-the-notebook"
        else:
            logger.important("source code is already saved")
    else:
        _source_code_artifact = ln.Artifact(
            _source_code_artifact_path,
            description=f"Source of transform {transform.uid}",
            version=transform.version,
            is_new_version_of=prev_source,
            visibility=0,  # hidden file
            run=False,
        )
        _source_code_artifact.save(upload=True, print_progress=False)
        transform._source_code_artifact = _source_code_artifact
        logger.debug(
            f"saved transform._source_code_artifact: {transform._source_code_artifact}"
        )

    # track environment
    env_path = ln_setup.settings.storage.cache_dir / f"run_env_pip_{run.uid}.txt"
    if env_path.exists():
        overwrite_env = True
        if run.environment_id is not None and from_cli:
            logger.important("run.environment is already saved")
            overwrite_env = False
        if overwrite_env:
            hash, _ = hash_file(env_path)
            artifact = ln.Artifact.filter(hash=hash, visibility=0).one_or_none()
            new_env_artifact = artifact is None
            if new_env_artifact:
                artifact = ln.Artifact(
                    env_path,
                    description="requirements.txt",
                    visibility=0,
                    run=False,
                )
                artifact.save(upload=True, print_progress=False)
            run.environment = artifact
            if new_env_artifact:
                logger.debug(f"saved run.environment: {run.environment}")

    # set finished_at
    if finished_at:
        run.finished_at = datetime.now(timezone.utc)

    # track report and set is_consecutive
    if not is_notebook:
        run.is_consecutive = True
        run.save()
    else:
        if run.report_id is not None:
            hash, _ = hash_file(report_path)  # ignore hash_type for now
            if hash != run.report.hash:
                if os.getenv("LAMIN_TESTING") is None:
                    # in test, auto-confirm overwrite
                    response = input(
                        f"You are about to replace (overwrite) an existing run report (hash '{run.report.hash}'). Proceed? (y/n)"
                    )
                else:
                    response = "y"
                if response == "y":
                    run.report.replace(report_path)
                    run.report.save(upload=True)
                else:
                    logger.important("keeping old report")
            else:
                logger.important("report is already saved")
        else:
            report_file = ln.Artifact(
                report_path,
                description=f"Report of run {run.uid}",
                is_new_version_of=prev_report,
                visibility=0,  # hidden file
                run=False,
            )
            report_file.save(upload=True, print_progress=False)
            run.report = report_file
        run.is_consecutive = is_consecutive
        run.save()
        logger.debug(
            f"saved transform.latest_run.report: {transform.latest_run.report}"
        )
    transform.save()

    # finalize
    if ln_setup.settings.instance.is_on_hub:
        identifier = ln_setup.settings.instance.slug
        logger.important(
            f"go to: https://lamin.ai/{identifier}/transform/{transform.uid}"
        )
        if not from_cli:
            thing, name = (
                ("notebook", "notebook.ipynb")
                if is_notebook
                else ("script", "script.py")
            )
            logger.important(
                f"if you want to update your {thing} without re-running it, use `lamin save {name}`"
            )
    # because run & transform changed, update the global context
    context.run = run
    context.run.transform = transform
    return None
