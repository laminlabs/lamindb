from __future__ import annotations

import os
import shutil
import subprocess
from datetime import datetime, timezone
from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
from lamin_utils import logger
from lnschema_core.types import TransformType

from .core._run_context import is_run_from_ipython, run_context

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


def finish():
    """Mark a tracked run as finished.

    If run in a notebook, it saves the run report & source code to your default storage location.
    """
    if run_context.path is None:
        raise TrackNotCalled("Please run `ln.track()` before `ln.finish()`")
    if is_run_from_ipython:  # notebooks
        if (
            get_seconds_since_modified(run_context.path) > 3
            and os.getenv("LAMIN_TESTING") is None
        ):
            raise NotebookNotSaved(
                "Please save the notebook in your editor right before running `ln.finish()`"
            )
        save_run_context_core(
            run=run_context.run,
            transform=run_context.transform,
            filepath=run_context.path,
            finished_at=True,
        )
    else:  # scripts
        # save_run_context_core was already called during ln.track()
        run_context.run.finished_at = datetime.now(timezone.utc)  # update run time
        run_context.run.save()


def save_run_context_core(
    *,
    run: Run,
    transform: Transform,
    filepath: Path,
    transform_family: QuerySet | None = None,
    finished_at: bool = False,
) -> str | None:
    import lamindb as ln

    ln.settings.verbosity = "success"

    # for scripts, things are easy
    is_consecutive = True
    source_code_path = filepath
    # for notebooks, we need more work
    if transform.type == TransformType.notebook:
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
        is_consecutive = check_consecutiveness(notebook_content)
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
            "jupyter nbconvert --to html"
            f" '{filepath.as_posix()}' --Application.log_level=40",
            shell=True,
            check=True,
        )
        # move the temporary file into the cache dir in case it's accidentally
        # in an existing storage location -> we want to move associated
        # artifacts into default storage and not register them in an existing
        # location
        filepath_html_orig = filepath.with_suffix(".html")  # current location
        filepath_html = ln_setup.settings.storage.cache_dir / filepath_html_orig.name
        # don't use Path.rename here because of cross-device link error
        # https://laminlabs.slack.com/archives/C04A0RMA0SC/p1710259102686969
        shutil.move(
            filepath_html_orig,  # type: ignore
            filepath_html,
        )
        # strip the output from the notebook to create the source code file
        # first, copy the notebook file to a temporary file in the cache
        source_code_path = ln_setup.settings.storage.cache_dir / filepath.name
        shutil.copy2(filepath, source_code_path)  # copy
        subprocess.run(
            f"nbstripout '{source_code_path}' --extra-keys='metadata.version metadata.kernelspec metadata.language_info metadata.pygments_lexer metadata.name metadata.file_extension'",
            shell=True,
            check=True,
        )
    # find initial versions of source codes and html reports
    prev_report = None
    prev_source = None
    if transform_family is None:
        transform_family = transform.versions
    if len(transform_family) > 0:
        for prev_transform in transform_family.order_by("-created_at"):
            if prev_transform.latest_report_id is not None:
                prev_report = prev_transform.latest_report
            if prev_transform.source_code_id is not None:
                prev_source = prev_transform.source_code
    ln.settings.silence_file_run_transform_warning = True
    # register the source code
    if transform.source_code is not None:
        # check if the hash of the notebook source code matches
        check_source_code = ln.Artifact(source_code_path, key="dummy")
        if check_source_code._state.adding:
            if os.getenv("LAMIN_TESTING") is None:
                # in test, auto-confirm overwrite
                response = input(
                    f"You are about to overwrite existing source code (hash {transform.source_code.hash}) for transform version"
                    f" '{transform.version}'. Proceed? (y/n)"
                )
            else:
                response = "y"
            if response == "y":
                transform.source_code.replace(source_code_path)
                transform.source_code.save(upload=True)
            else:
                logger.warning("Please re-run `ln.track()` to make a new version")
                return "rerun-the-notebook"
    else:
        source_code = ln.Artifact(
            source_code_path,
            description=f"Source of transform {transform.uid}",
            version=transform.version,
            is_new_version_of=prev_source,
            visibility=0,  # hidden file
            run=False,
        )
        source_code.save(upload=True)
        transform.source_code = source_code
        logger.success(f"saved transform.source_code: {transform.source_code}")
    # track environment
    filepath_env = ln_setup.settings.storage.cache_dir / f"run_env_pip_{run.uid}.txt"
    if filepath_env.exists():
        artifact = ln.Artifact(
            filepath_env,
            description="requirements.txt",
            visibility=0,
            run=False,
        )
        if artifact._state.adding:
            artifact.save(upload=True)
        run.environment = artifact
        logger.success(f"saved run.environment: {run.environment}")
    # save report file
    if not transform.type == TransformType.notebook:
        run.save()
    else:
        if run.report_id is not None:
            logger.warning(
                "there is already an existing report for this run, replacing it"
            )
            run.report.replace(filepath_html)
            run.report.save(upload=True)
        else:
            report_file = ln.Artifact(
                filepath_html,
                description=f"Report of run {run.uid}",
                is_new_version_of=prev_report,
                visibility=0,  # hidden file
                run=False,
            )
            report_file.save(upload=True)
            run.report = report_file
        run.is_consecutive = is_consecutive
        if finished_at:
            run.finished_at = datetime.now(timezone.utc)
        run.save()
        transform.latest_report = run.report
    transform.save()
    if transform.type == TransformType.notebook:
        logger.success(f"saved transform.latest_report: {transform.latest_report}")
    if ln_setup.settings.instance.is_remote:
        identifier = ln_setup.settings.instance.slug
        logger.success(
            f"go to: https://lamin.ai/{identifier}/transform/{transform.uid}"
        )
    # because run & transform changed, update the global run_context
    run_context.run = run
    run_context.transform = transform
    return None
