from __future__ import annotations

import re
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
        # strip resaved finish error if present
        # this is normally the last cell
        if cell["cell_type"] == "code" and ".finish(" in cell["source"]:
            for output in cell["outputs"]:
                if output.get("ename", None) == "NotebookNotSaved":
                    cell["outputs"] = []
                    break
    return None


def notebook_to_report(notebook_path: Path, output_path: Path) -> None:
    import nbformat
    import traitlets.config as config
    from nbconvert import HTMLExporter

    with open(notebook_path, encoding="utf-8") as f:
        notebook = nbformat.read(f, as_version=4)
    prepare_notebook(notebook, strip_title=True)
    notebook.metadata.clear()  # strip notebook metadata
    # if we were to export as ipynb, the following two lines would do it
    # with open(output_path, "w", encoding="utf-8") as f:
    #     nbformat.write(notebook, f)
    # instead we need all this code
    c = config.Config()
    c.HTMLExporter.preprocessors = []
    c.HTMLExporter.exclude_input_prompt = True
    c.HTMLExporter.exclude_output_prompt = True
    c.HTMLExporter.anchor_link_text = " "
    html_exporter = HTMLExporter(config=c)
    html, _ = html_exporter.from_notebook_node(notebook)
    output_path.write_text(html, encoding="utf-8")


def notebook_to_script(
    transform: Transform, notebook_path: Path, script_path: Path
) -> None:
    import jupytext

    notebook = jupytext.read(notebook_path)
    py_content = jupytext.writes(notebook, fmt="py:percent")
    # remove global metadata header
    py_content = re.sub(r"^# ---\n.*?# ---\n\n", "", py_content, flags=re.DOTALL)
    # replace title
    py_content = py_content.replace(f"# # {transform.name}", "# # transform.name")
    script_path.write_text(py_content)


def save_context_core(
    *,
    run: Run,
    transform: Transform,
    filepath: Path,
    finished_at: bool = False,
    ignore_non_consecutive: bool | None = None,
    from_cli: bool = False,
) -> str | None:
    from lnschema_core.models import (
        format_field_value,  # needs to come after lamindb was imported because of CLI use
    )

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
        if not ignore_non_consecutive:  # ignore_non_consecutive is None or False
            is_consecutive = check_consecutiveness(
                notebook_content, calling_statement=".finish("
            )
            if not is_consecutive:
                response = "n"  # ignore_non_consecutive == False
                if ignore_non_consecutive is None:
                    response = input(
                        "   Do you still want to proceed with finishing? (y/n) "
                    )
                if response != "y":
                    return "aborted-non-consecutive"
        # write the report
        report_path = ln_setup.settings.cache_dir / filepath.name.replace(
            ".ipynb", ".html"
        )
        notebook_to_report(filepath, report_path)
        # write the source code
        source_code_path = ln_setup.settings.cache_dir / filepath.name.replace(
            ".ipynb", ".py"
        )
        notebook_to_script(transform, filepath, source_code_path)
    ln.settings.creation.artifact_silence_missing_run_warning = True
    # track source code
    hash, _ = hash_file(source_code_path)  # ignore hash_type for now
    if (
        transform._source_code_artifact_id is not None
        or transform.source_code is not None  # equivalent to transform.hash is not None
    ):
        # check if the hash of the transform source code matches
        # (for scripts, we already run the same logic in track() - we can deduplicate the call at some point)
        ref_hash = (
            transform.hash
            if transform.hash is not None
            else transform._source_code_artifact.hash
        )
        if hash != ref_hash:
            response = input(
                f"You are about to overwrite existing source code (hash '{ref_hash}') for Transform('{transform.uid}')."
                f" Proceed? (y/n)"
            )
            if response == "y":
                transform.source_code = source_code_path.read_text()
                transform.hash = hash
            else:
                logger.warning("Please re-run `ln.track()` to make a new version")
                return "rerun-the-notebook"
        else:
            logger.important("source code is already saved")
    else:
        transform.source_code = source_code_path.read_text()
        transform.hash = hash

    # track environment
    env_path = ln_setup.settings.cache_dir / f"run_env_pip_{run.uid}.txt"
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
                response = input(
                    f"You are about to overwrite an existing report (hash '{run.report.hash}') for Run('{run.uid}'). Proceed? (y/n)"
                )
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
    if not from_cli:
        run_time = run.finished_at - run.started_at
        logger.important(
            f"finished Run('{run.uid[:8]}') after {run_time} at {format_field_value(run.finished_at)}"
        )
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
