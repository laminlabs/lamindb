import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import lamindb_setup as ln_setup
from lamin_utils import logger
from lnschema_core import Run, Transform
from lnschema_core.types import TransformType

from ._query_set import QuerySet
from .core._run_context import is_run_from_ipython, run_context


class CallFinishInLastCell(SystemExit):
    pass


def finish(i_saved_the_notebook: bool = False):
    """Mark a tracked run as finished.

    When run in notebooks, save the run report to your default storage location.

    Args:
        i_saved_the_notebook: Indicate that you saved the notebook in your
            editor (JupyterLab, VSCode, etc.).
    """
    if is_run_from_ipython:
        # notebooks
        from nbproject.dev import read_notebook
        from nbproject.dev._check_last_cell import check_last_cell

        if not i_saved_the_notebook and not ln_setup._TESTING:
            logger.error(
                "Please save the notebook, pass `i_saved_the_notebook=True`, and re-run this cell."
            )
            return None
        notebook_content = read_notebook(run_context.path)  # type: ignore
        if not check_last_cell(notebook_content, "i_saved_the_notebook"):
            raise CallFinishInLastCell(
                "Can only finish() from the last code cell of the notebook."
            )
        save_run_context_core(
            run=run_context.run,
            transform=run_context.transform,
            filepath=run_context.path,
            finished_at=True,
            notebook_content=notebook_content,
        )
    else:
        # scripts
        run_context.run.finished_at = datetime.now(timezone.utc)  # update run time
        run_context.run.save()


# do not type because we need to be aware of lnschema_core import order
def save_run_context_core(
    *,
    run: Run,
    transform: Transform,
    filepath: Path,
    transform_family: Optional[QuerySet] = None,
    is_consecutive: bool = True,
    finished_at: bool = False,
    notebook_content=None,  # nbproject.Notebook
) -> Optional[str]:
    import lamindb as ln

    ln.settings.verbosity = "success"

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
        if notebook_content is None:
            notebook_content = read_notebook(filepath)  # type: ignore
        is_consecutive = check_consecutiveness(notebook_content)
        if not is_consecutive:
            if os.getenv("LAMIN_TESTING") is None:
                decide = input(
                    "   Do you still want to proceed with publishing? (y/n) "
                )
            else:
                decide = "n"
            if decide != "y":
                logger.error("Aborted (non-consecutive)!")
                return "aborted-non-consecutive"

        # convert the notebook file to html
        # log_level is set to 40 to silence the nbconvert logging
        result = subprocess.run(
            "jupyter nbconvert --to html"
            f" {filepath.as_posix()} --Application.log_level=40",
            shell=True,
        )
        # move the temporary file into the cache dir in case it's accidentally
        # in an existing storage location -> we want to move associated
        # artifacts into default storage and not register them in an existing
        # location
        filepath_html = filepath.with_suffix(".html")  # current location
        shutil.move(
            filepath_html,  # type: ignore
            ln_setup.settings.storage.cache_dir / filepath_html.name,
        )  # move; don't use Path.rename here because of cross-device link error
        # see https://laminlabs.slack.com/archives/C04A0RMA0SC/p1710259102686969
        filepath_html = (
            ln_setup.settings.storage.cache_dir / filepath_html.name
        )  # adjust location
        assert result.returncode == 0
        # copy the notebook file to a temporary file
        source_code_path = ln_setup.settings.storage.cache_dir / filepath.name
        shutil.copy2(filepath, source_code_path)  # copy
        result = subprocess.run(f"nbstripout {source_code_path}", shell=True)
        assert result.returncode == 0
    else:
        source_code_path = filepath
    # find initial versions of source codes and html reports
    initial_report = None
    initial_source = None
    if transform_family is None:
        transform_family = transform.versions
    if len(transform_family) > 0:
        for prev_transform in transform_family.order_by("-created_at"):
            # check for id to avoid query
            if prev_transform.latest_report_id is not None:
                # any previous latest report of this transform is OK!
                initial_report = prev_transform.latest_report
            if prev_transform.source_code_id is not None:
                # any previous source code id is OK!
                initial_source = prev_transform.source_code
    ln.settings.silence_file_run_transform_warning = True
    # register the source code
    if transform.source_code is not None:
        # check if the hash of the notebook source code matches
        check_source_code = ln.Artifact(source_code_path, key="dummy")
        if check_source_code._state.adding:
            if os.getenv("LAMIN_TESTING") is None:
                # in test, auto-confirm overwrite
                response = input(
                    "You try to save a new notebook source code with the same version"
                    f" '{transform.version}'; do you want to replace the content of the"
                    f" existing source code {transform.source_code}? (y/n)"
                )
            else:
                response = "y"
            if response == "y":
                transform.source_code.replace(source_code_path)
                transform.source_code.save()
            else:
                logger.warning(
                    "Please create a new version of the notebook via `lamin track"
                    " <filepath>` and re-run the notebook"
                )
                return "rerun-the-notebook"
    else:
        source_code = ln.Artifact(
            source_code_path,
            description=f"Source of transform {transform.uid}",
            version=transform.version,
            is_new_version_of=initial_source,
            visibility=0,  # hidden file
            run=False,
        )
        source_code.save()
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
            artifact.save()
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
            run.report.save()
        else:
            report_file = ln.Artifact(
                filepath_html,
                description=f"Report of run {run.uid}",
                is_new_version_of=initial_report,
                visibility=0,  # hidden file
                run=False,
            )
            report_file.save()
            run.report = report_file
        run.is_consecutive = is_consecutive
        if finished_at:
            run.finished_at = datetime.now(timezone.utc)
        run.save()
        transform.latest_report = run.report
    transform.save()
    if transform.type == TransformType.notebook:
        logger.success(f"saved transform.latest_report: {transform.latest_report}")
    identifier = ln_setup.settings.instance.slug
    logger.success(f"go to: https://lamin.ai/{identifier}/transform/{transform.uid}")
    # because run & transform changed, update the global run_context
    run_context.run = run
    run_context.transform = transform
    return None
