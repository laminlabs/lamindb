from __future__ import annotations

import os
import re
import shutil
from datetime import datetime, timezone
from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
from lamin_utils import logger
from lamindb_setup.core.hashing import hash_file

if TYPE_CHECKING:
    from pathlib import Path

    from lnschema_core import Run, Transform

    from ._query_set import QuerySet


# this is from the get_title function in nbproject
# should be moved into lamindb sooner or later
def prepare_notebook(
    nb,
    strip_title: bool = False,
) -> str | None:
    """Strip title from the notebook if requested."""
    title_found = False
    for cell in nb.cells:
        cell.metadata.clear()  # strip cell metadata
        if not title_found and cell["cell_type"] == "markdown":
            lines = cell["source"].split("\n")
            for i, line in enumerate(lines):
                if line.startswith("# "):
                    line.lstrip("#").strip(" .").strip()
                    title_found = True
                    if strip_title:
                        lines.pop(i)
                        cell["source"] = "\n".join(lines)
    return None


def notebook_to_ipynb(notebook_path: Path, output_path: Path) -> None:
    import nbformat

    with open(notebook_path, encoding="utf-8") as f:
        notebook = nbformat.read(f, as_version=4)
    prepare_notebook(notebook, strip_title=True)
    notebook.metadata.clear()  # strip notebook metadata
    with open(output_path, "w", encoding="utf-8") as f:
        nbformat.write(notebook, f)


def notebook_to_script(
    transform: Transform, notebook_path: Path, script_path: Path
) -> None:
    import jupytext

    notebook = jupytext.read(notebook_path)
    py_content = jupytext.writes(notebook, fmt="py:percent")
    # remove global metadata header
    py_content = re.sub(r"^# ---\n.*?# ---\n\n", "", py_content, flags=re.DOTALL)
    # replace title
    py_content = py_content.replace(f"# {transform.name}", "# Transform.name")
    script_path.write_text(py_content)


def script_to_notebook(transform: Transform, notebook_path: Path) -> None:
    import jupytext

    notebook = jupytext.reads(transform.source_code, fmt="py:percent")
    prepare_notebook(notebook)
    jupytext.write(notebook, notebook_path)


def save_context_core(
    *,
    run: Run,
    transform: Transform,
    filepath: Path,
    finished_at: bool = False,
    from_cli: bool = False,
) -> str | None:
    import lamindb as ln

    from .core._context import context, is_run_from_ipython

    ln.settings.verbosity = "success"

    # for scripts, things are easy
    is_consecutive = True
    is_notebook = transform.type == "notebook"
    source_code_path = filepath
    # for notebooks, we need more work
    if is_notebook:
        try:
            import jupytext
            from nbproject.dev import (
                check_consecutiveness,
                read_notebook,
            )
        except ImportError:
            logger.error("install nbproject & jupytext: pip install nbproject jupytext")
            return None
        notebook_content = read_notebook(filepath)  # type: ignore
        is_consecutive = check_consecutiveness(
            notebook_content, calling_statement=".finish()"
        )
        if not is_consecutive:
            msg = "   Do you still want to proceed with finishing? (y/n) "
            if os.getenv("LAMIN_TESTING") is None:
                response = input(msg)
            else:
                response = "n"
            if response != "y":
                return "aborted-non-consecutive"
        # write the report
        report_path = ln_setup.settings.storage.cache_dir / filepath.name
        notebook_to_ipynb(filepath, report_path)
        # write the source code
        source_code_path = ln_setup.settings.storage.cache_dir / filepath.name
        notebook_to_script(transform, filepath, source_code_path)
    ln.settings.creation.artifact_silence_missing_run_warning = True
    # track source code
    hash, _ = hash_file(source_code_path)  # ignore hash_type for now
    if (
        transform._source_code_artifact_id is not None
        or transform.source_code is not None
    ):
        # check if the hash of the transform source code matches
        # (for scripts, we already run the same logic in track() - we can deduplicate the call at some point)
        if transform.hash is not None:
            condition = hash != transform.hash
        else:
            condition = hash != transform._source_code_artifact.hash
        if condition:
            if os.getenv("LAMIN_TESTING") is None:
                # in test, auto-confirm overwrite
                response = input(
                    f"You are about to replace (overwrite) existing source code (hash '{transform._source_code_artifact.hash}') for transform version"
                    f" '{transform.version}'. Proceed? (y/n)"
                )
            else:
                response = "y"
            if response == "y":
                transform.source_code = source_code_path.read_text()
                transform.hash = hash
            else:
                logger.warning(
                    "Please re-run `ln.context.track()` to make a new version"
                )
                return "rerun-the-notebook"
        else:
            logger.important("source code is already saved")
    else:
        transform.source_code = source_code_path.read_text()
        transform.hash = hash

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
    context._run = run
    context._transform = transform
    return None
