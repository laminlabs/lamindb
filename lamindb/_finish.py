from __future__ import annotations

import builtins
import re
from datetime import datetime, timezone
from time import sleep
from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
from lamin_utils import logger
from lamin_utils._logger import LEVEL_TO_COLORS, LEVEL_TO_ICONS, RESET_COLOR
from lamindb_setup.core.hashing import hash_file

from lamindb.models import Artifact, Run, Transform

is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)

if TYPE_CHECKING:
    from pathlib import Path


def get_save_notebook_message() -> str:
    # do not add bold() or any other complicated characters as then we can't match this
    # easily anymore in an html to strip it out
    return f"please hit {get_shortcut()} to save the notebook in your editor"


def get_save_notebook_message_retry() -> str:
    return f"{get_save_notebook_message()} and re-run finish()"


# this code was originally in nbproject by the same authors
def check_consecutiveness(
    nb, calling_statement: str = None, silent_success: bool = True
) -> bool:
    """Check whether code cells have been executed consecutively.

    Needs to be called in the last code cell of a notebook.
    Otherwise raises `RuntimeError`.

    Returns cell transitions that violate execution at increments of 1 as a list
    of tuples.

    Args:
        nb: Notebook content.
        calling_statement: The statement that calls this function.
    """
    cells = nb.cells

    violations = []
    prev = 0

    ccount = 0  # need to initialize because notebook might note have code cells
    # and below, we check if ccount is None
    for cell in cells:
        cell_source = "".join(cell["source"])
        if cell["cell_type"] != "code" or cell_source == "":
            continue

        if calling_statement is not None and calling_statement in cell_source:
            continue

        ccount = cell["execution_count"]
        if ccount is None or prev is None or ccount - prev != 1:
            violations.append((prev, ccount))

        prev = ccount

    # ignore the very last code cell of the notebook
    # `check_consecutiveness` is being run during publish if `last_cell`` is True
    # hence, that cell has ccount is None
    if ccount is None:
        violations.pop()

    any_violations = len(violations) > 0
    if any_violations:
        logger.warning(f"cells {violations} were not run consecutively")
    elif not silent_success:
        logger.success("cell execution numbers increase consecutively")

    return not any_violations


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
        artifact = Artifact(  # type: ignore
            logs_path,
            description=f"log streams of run {run.uid}",
            kind="__lamindb_run__",
            run=False,
        )
        artifact.save(upload=True, print_progress=False)
        run.report = artifact
        if save_run:  # defaults to false because is slow
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
        # strip logging message about saving notebook in editor
        # this is normally the last cell
        if cell["cell_type"] == "code" and ".finish(" in cell["source"]:
            for output in cell["outputs"]:
                if "to save the notebook in your editor" in output.get("text", ""):
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


def notebook_to_script(  # type: ignore
    title: str, notebook_path: Path, script_path: Path | None = None
) -> None | str:
    import jupytext

    notebook = jupytext.read(notebook_path)
    py_content = jupytext.writes(notebook, fmt="py:percent")
    # remove global metadata header
    py_content = re.sub(r"^# ---\n.*?# ---\n\n", "", py_content, flags=re.DOTALL)
    # replace title
    py_content = py_content.replace(f"# # {title}", "#")
    if script_path is None:
        return py_content
    else:
        script_path.write_text(py_content)


def clean_r_notebook_html(file_path: Path) -> tuple[str | None, Path]:
    import re

    cleaned_content = file_path.read_text()
    # remove title from content
    pattern_title = r"<title>(.*?)</title>"
    title_match = re.search(pattern_title, cleaned_content)
    title_text = None
    if title_match:
        title_text = title_match.group(1)
        pattern_h1 = f"<h1[^>]*>{re.escape(title_text)}</h1>"
        cleaned_content = re.sub(pattern_title, "", cleaned_content)
        cleaned_content = re.sub(pattern_h1, "", cleaned_content)
    # remove error message from content
    if "to save the notebook in your editor" in cleaned_content:
        orig_error_message = f"! {get_save_notebook_message_retry()}"
        # coming up with the regex for this is a bit tricky due to all the
        # escape characters we'd need to insert into the message; hence,
        # we do this with a replace() instead
        cleaned_content = cleaned_content.replace(orig_error_message, "")
        if "to save the notebook in your editor" in cleaned_content:
            orig_error_message = orig_error_message.replace(
                " finish()", "\nfinish()"
            )  # RStudio might insert a newline
            cleaned_content = cleaned_content.replace(orig_error_message, "")
    cleaned_path = file_path.parent / (f"{file_path.stem}.cleaned{file_path.suffix}")
    cleaned_path.write_text(cleaned_content)
    return title_text, cleaned_path


def check_filepath_recently_saved(filepath: Path, is_finish_retry: bool) -> bool:
    # the recently_saved_time needs to be very low for the first check
    # because an accidental save (e.g. via auto-save) might otherwise lead
    # to upload of an outdated notebook
    # also see implementation for R notebooks below
    offset_saved_time = 0.3 if not is_finish_retry else 20
    for retry in range(30):
        recently_saved_time = offset_saved_time + retry  # sleep time is 1 sec
        if get_seconds_since_modified(filepath) > recently_saved_time:
            if retry == 0:
                prefix = f"{LEVEL_TO_COLORS[20]}{LEVEL_TO_ICONS[20]}{RESET_COLOR}"
                print(f"{prefix} {get_save_notebook_message()}", end=" ")
            elif retry == 9:
                print(".", end="\n")
            elif retry == 4:
                print(". still waiting ", end="")
            else:
                print(".", end="")
            sleep(1)
        else:
            if retry > 0:
                prefix = f"{LEVEL_TO_COLORS[25]}{LEVEL_TO_ICONS[25]}{RESET_COLOR}"
                print(f" {prefix}")
            # filepath was recently saved, return True
            return True
    # if we arrive here, no save event occured, return False
    return False


def save_context_core(
    *,
    run: Run | None,
    transform: Transform,
    filepath: Path,
    finished_at: bool = False,
    ignore_non_consecutive: bool | None = None,
    from_cli: bool = False,
    is_retry: bool = False,
    notebook_runner: str | None = None,
) -> str | None:
    import lamindb as ln
    from lamindb.models import (
        format_field_value,  # needs to come after lamindb was imported because of CLI use
    )

    ln.settings.verbosity = "success"

    # for scripts, things are easy
    is_consecutive = True
    is_ipynb = filepath.suffix == ".ipynb"
    is_r_notebook = filepath.suffix in {".qmd", ".Rmd"}
    source_code_path = filepath
    report_path: Path | None = None
    save_source_code_and_report = True
    if (
        is_run_from_ipython and notebook_runner != "nbconvert"
    ):  # python notebooks in interactive session
        import nbproject

        # it might be that the user modifies the title just before ln.finish()
        if (nbproject_title := nbproject.meta.live.title) != transform.description:
            transform.description = nbproject_title
            transform.save()
        if not ln_setup._TESTING:
            save_source_code_and_report = check_filepath_recently_saved(
                filepath, is_retry
            )
            if not save_source_code_and_report and not is_retry:
                logger.warning(get_save_notebook_message_retry())
                return "retry"
            elif not save_source_code_and_report:
                logger.warning(
                    "the notebook on disk wasn't saved within the last 10 sec"
                )
    if is_ipynb:  # could be from CLI outside interactive session
        try:
            import jupytext  # noqa: F401
            from nbproject.dev import (
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
                if ignore_non_consecutive is None:  # only print warning
                    response = "y"  # we already printed the warning
                else:  # ask user to confirm
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
        notebook_to_script(transform.description, filepath, source_code_path)
    elif is_r_notebook:
        if filepath.with_suffix(".nb.html").exists():
            report_path = filepath.with_suffix(".nb.html")
        elif filepath.with_suffix(".html").exists():
            report_path = filepath.with_suffix(".html")
        else:
            logger.warning(
                f"no html report found; to attach one, create an .html export for your {filepath.suffix} file and then run: lamin save {filepath}"
            )
    if report_path is not None and is_r_notebook and not from_cli:  # R notebooks
        # see comment above in check_filepath_recently_saved
        recently_saved_time = 0.3 if not is_retry else 20
        if get_seconds_since_modified(report_path) > recently_saved_time:
            # the automated retry solution of Jupyter notebooks does not work in RStudio because the execution of the notebook cell
            # seems to block the event loop of the frontend
            if not is_retry:
                logger.warning(get_save_notebook_message_retry())
                return "retry"
            else:
                logger.warning(
                    "the notebook on disk hasn't been saved within the last 20 sec"
                )
            save_source_code_and_report = False
    ln.settings.creation.artifact_silence_missing_run_warning = True
    # save source code
    if save_source_code_and_report:
        transform_hash, _ = hash_file(source_code_path)  # ignore hash_type for now
        if transform.hash is not None:
            # check if the hash of the transform source code matches
            # (for scripts, we already run the same logic in track() - we can deduplicate the call at some point)
            if transform_hash != transform.hash:
                response = input(
                    f"You are about to overwrite existing source code (hash '{transform.hash}') for Transform('{transform.uid}')."
                    f" Proceed? (y/n)"
                )
                if response == "y":
                    transform.source_code = source_code_path.read_text()
                    transform.hash = transform_hash
                else:
                    logger.warning("Please re-run `ln.track()` to make a new version")
                    return "rerun-the-notebook"
            else:
                logger.debug("source code is already saved")
        else:
            transform.source_code = source_code_path.read_text()
            transform.hash = transform_hash

    # track run environment
    if run is not None:
        env_path = ln_setup.settings.cache_dir / f"run_env_pip_{run.uid}.txt"
        if env_path.exists():
            overwrite_env = True
            if run.environment_id is not None and from_cli:
                logger.important("run.environment is already saved, ignoring")
                overwrite_env = False
            if overwrite_env:
                env_hash, _ = hash_file(env_path)
                artifact = ln.Artifact.objects.filter(hash=env_hash).one_or_none()
                new_env_artifact = artifact is None
                if new_env_artifact:
                    artifact = ln.Artifact(  # type: ignore
                        env_path,
                        description="requirements.txt",
                        kind="__lamindb_run__",
                        run=False,
                    )
                    artifact.save(upload=True, print_progress=False)
                run.environment = artifact
                if new_env_artifact:
                    logger.debug(f"saved run.environment: {run.environment}")

    # set finished_at
    if finished_at and run is not None:
        if not from_cli:
            update_finished_at = True
        else:
            update_finished_at = run.finished_at is None
        if update_finished_at:
            run.finished_at = datetime.now(timezone.utc)

    # track logs
    if run is not None and not from_cli and not is_ipynb and not is_r_notebook:
        save_run_logs(run)

    # track report and set is_consecutive
    if save_source_code_and_report:
        if run is not None:
            # do not save a run report if executing through nbconvert
            if report_path is not None and notebook_runner != "nbconvert":
                if is_r_notebook:
                    title_text, report_path = clean_r_notebook_html(report_path)
                    if title_text is not None:
                        transform.description = title_text
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
                    report_file = ln.Artifact(  # type: ignore
                        report_path,
                        description=f"Report of run {run.uid}",
                        kind="__lamindb_run__",  # hidden file
                        run=False,
                    )
                    report_file.save(upload=True, print_progress=False)
                    run.report = report_file
                if is_r_notebook:
                    # this is the "cleaned" report
                    report_path.unlink()
                logger.debug(
                    f"saved transform.latest_run.report: {transform.latest_run.report}"
                )
            run._is_consecutive = is_consecutive
        if report_path is not None and notebook_runner == "nbconvert":
            logger.important(f"to save the notebook html, run: lamin save {filepath}")

    # save both run & transform records if we arrive here
    if run is not None:
        run.save()
    transform_id_prior_to_save = transform.id
    transform.save()  # this in-place updates the state of transform upon hash collision
    if transform.id != transform_id_prior_to_save:
        # the hash existed and we're actually back to the previous version
        # hence, this was in fact a run of the previous transform rather than of
        # the new transform
        # this can happen in interactively executed notebooks with a pro-active version bump in case it turns out that the user didn't make a change to the notebook
        run.transform = transform
        run.save()
        ln.Transform.get(transform_id_prior_to_save).delete()

    # finalize
    if not from_cli and run is not None:
        run_time = run.finished_at - run.started_at
        days = run_time.days
        seconds = run_time.seconds
        hours = seconds // 3600
        minutes = (seconds % 3600) // 60
        secs = seconds % 60
        formatted_run_time = (
            f"{days}d"
            if days != 0
            else "" + f"{hours}h"
            if hours != 0
            else "" + f"{minutes}m"
            if minutes != 0
            else "" + f"{secs}s"
        )

        logger.important(
            f"finished Run('{run.uid[:8]}') after {formatted_run_time} at {format_field_value(run.finished_at)}"
        )
    if ln_setup.settings.instance.is_on_hub:
        instance_slug = ln_setup.settings.instance.slug
        if save_source_code_and_report:
            logger.important(
                f"go to: https://lamin.ai/{instance_slug}/transform/{transform.uid}"
            )
        if not from_cli and save_source_code_and_report:
            thing = "notebook" if (is_ipynb or is_r_notebook) else "script"
            logger.important(
                f"to update your {thing} from the CLI, run: lamin save {filepath}"
            )
    if not save_source_code_and_report:
        logger.warning(
            f"did *not* save source code and report -- to do so, run: lamin save {filepath}"
        )
    return None
