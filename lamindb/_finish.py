from __future__ import annotations

import re
from datetime import datetime, timezone
from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
from lamin_utils import logger
from lamindb_setup.core.hashing import hash_file

from lamindb.core.exceptions import NotebookNotSaved
from lamindb.models import Artifact, Run, Transform

if TYPE_CHECKING:
    from pathlib import Path


def get_save_notebook_message() -> str:
    return f"Please save the notebook in your editor (shortcut `{get_shortcut()}`) within 2 sec before calling `finish()`"


def get_shortcut() -> str:
    import platform

    return "CMD + s" if platform.system() == "Darwin" else "CTRL + s"


def get_seconds_since_modified(filepath) -> float:
    return datetime.now().timestamp() - filepath.stat().st_mtime


def save_run_logs(run: Run, save_run: bool = False) -> None:
    logs_path = ln_setup.settings.cache_dir / f"run_logs_{run.uid}.txt"
    if logs_path.exists():
        if run.report is not None:
            logger.important("overwriting run.report")
        artifact = Artifact(
            logs_path,
            description=f"log streams of run {run.uid}",
            visibility=0,
            run=False,
        )
        artifact.save(upload=True, print_progress=False)
        run.report = artifact
        if save_run:  # defaults to fast because is slow
            run.save()


# this is from the get_title function in nbproject
# should be moved into lamindb sooner or later
def prepare_notebook(
    nb,
    strip_title: bool = False,
) -> str | None:
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


# removes NotebookNotSaved error message from notebook html
def clean_r_notebook_html(file_path: Path) -> tuple[str | None, Path]:
    import re

    cleaned_content = (
        file_path.read_text()
    )  # at this point cleaned_content is still raw
    pattern_title = r"<title>(.*?)</title>"
    title_match = re.search(pattern_title, cleaned_content)
    title_text = None
    if title_match:
        title_text = title_match.group(1)
        pattern_h1 = f"<h1[^>]*>{re.escape(title_text)}</h1>"
        cleaned_content = re.sub(pattern_title, "", cleaned_content)
        cleaned_content = re.sub(pattern_h1, "", cleaned_content)
    cleaned_content = cleaned_content.replace(
        f"NotebookNotSaved: {get_save_notebook_message()}", ""
    )
    cleaned_path = file_path.parent / (f"{file_path.stem}.cleaned{file_path.suffix}")
    cleaned_path.write_text(cleaned_content)
    return title_text, cleaned_path


def save_context_core(
    *,
    run: Run | None,
    transform: Transform,
    filepath: Path,
    finished_at: bool = False,
    ignore_non_consecutive: bool | None = None,
    from_cli: bool = False,
) -> str | None:
    import lamindb as ln
    from lamindb.models import (
        format_field_value,  # needs to come after lamindb was imported because of CLI use
    )

    from .core._context import context, is_run_from_ipython

    ln.settings.verbosity = "success"

    # for scripts, things are easy
    is_consecutive = True
    is_ipynb = filepath.suffix == ".ipynb"
    is_r_notebook = filepath.suffix in {".qmd", ".Rmd"}
    source_code_path = filepath
    report_path: Path | None = None
    # for notebooks, we need more work
    if is_ipynb:
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
    elif is_r_notebook:
        if filepath.with_suffix(".nb.html").exists():
            report_path = filepath.with_suffix(".nb.html")
        elif filepath.with_suffix(".html").exists():
            report_path = filepath.with_suffix(".html")
        else:
            logger.warning(
                f"no {filepath.with_suffix('.nb.html')} found, save your manually rendered .html report via the CLI: lamin save {filepath}"
            )
    if report_path is not None and not from_cli:
        if get_seconds_since_modified(report_path) > 2 and not ln_setup._TESTING:
            # this can happen when auto-knitting an html with RStudio
            raise NotebookNotSaved(get_save_notebook_message())
    ln.settings.creation.artifact_silence_missing_run_warning = True
    # track source code
    hash, _ = hash_file(source_code_path)  # ignore hash_type for now
    if (
        transform._source_code_artifact_id is not None
        or transform.hash is not None  # .hash is equivalent to .transform
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
            logger.debug("source code is already saved")
    else:
        transform.source_code = source_code_path.read_text()
        transform.hash = hash

    # track environment
    if run is not None:
        env_path = ln_setup.settings.cache_dir / f"run_env_pip_{run.uid}.txt"
        if env_path.exists():
            overwrite_env = True
            if run.environment_id is not None and from_cli:
                logger.important("run.environment is already saved, ignoring")
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
    if finished_at and run is not None:
        run.finished_at = datetime.now(timezone.utc)

    # track logs
    if run is not None and not from_cli and not is_ipynb and not is_r_notebook:
        save_run_logs(run)

    # track report and set is_consecutive
    if run is not None:
        if report_path is not None:
            if is_r_notebook:
                title_text, report_path = clean_r_notebook_html(report_path)
                if title_text is not None:
                    transform.name = title_text
            if run.report_id is not None:
                hash, _ = hash_file(report_path)  # ignore hash_type for now
                if hash != run.report.hash:
                    response = input(
                        f"You are about to overwrite an existing report (hash '{run.report.hash}') for Run('{run.uid}'). Proceed? (y/n)"
                    )
                    if response == "y":
                        run.report.replace(report_path)
                        run.report.save(upload=True, print_progress=False)
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
            logger.debug(
                f"saved transform.latest_run.report: {transform.latest_run.report}"
            )
        run.is_consecutive = is_consecutive

        # save both run & transform records if we arrive here
        run.save()
    transform.save()

    # finalize
    if not from_cli and run is not None:
        run_time = run.finished_at - run.started_at
        days = run_time.days
        seconds = run_time.seconds
        hours = seconds // 3600
        minutes = (seconds % 3600) // 60
        secs = seconds % 60
        formatted_run_time = f"{days}d {hours}h {minutes}m {secs}s"

        logger.important(
            f"finished Run('{run.uid[:8]}') after {formatted_run_time} at {format_field_value(run.finished_at)}"
        )
    if ln_setup.settings.instance.is_on_hub:
        identifier = ln_setup.settings.instance.slug
        logger.important(
            f"go to: https://lamin.ai/{identifier}/transform/{transform.uid}"
        )
        if not from_cli:
            thing = "notebook" if (is_ipynb or is_r_notebook) else "script"
            logger.important(
                f"if you want to update your {thing} without re-running it, use `lamin save {filepath}`"
            )
    # because run & transform changed, update the global context
    context._run = run
    context._transform = transform
    return None
