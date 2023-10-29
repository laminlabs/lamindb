import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

import lamindb_setup
from lamin_utils import colors, logger


# also see lamindb.dev._run_context.reinitialize_notebook for related code
def update_notebook_metadata(nb, notebook_path):
    from nbproject.dev import write_notebook
    from nbproject.dev._initialize import nbproject_id

    stem_id = nb.metadata["nbproject"]["id"]
    current_version = nb.metadata["nbproject"]["version"]
    logger.info(
        f"the notebook {notebook_path} is already tracked (stem_id='{stem_id}',"
        f" version: '{current_version}')"
    )
    updated = False
    # ask for generating new id
    if os.getenv("LAMIN_TESTING") is None:
        response = input("Do you want to generate a new id prefix? (y/n) ")
    else:
        response = "y"
    if response == "y":
        nb.metadata["nbproject"]["id"] = nbproject_id()
        updated = True
    else:
        response = input(
            f"The current version is '{current_version}' - do you want to set a new"
            " version? (y/n) "
        )
        if response == "y":
            new_version = input("Please type the version: ")
            nb.metadata["nbproject"]["version"] = new_version
            updated = True
    if updated:
        logger.save("updated notebook metadata")
        write_notebook(nb, notebook_path)


def track(notebook_path: str, pypackage: Optional[str] = None) -> None:
    try:
        from nbproject.dev import initialize_metadata, read_notebook, write_notebook
    except ImportError:
        logger.error("install nbproject: pip install nbproject")
        return None

    nb = read_notebook(notebook_path)
    if "nbproject" not in nb.metadata:
        if pypackage is not None:
            pypackage = [pp for pp in pypackage.split(",") if len(pp) > 0]  # type: ignore # noqa
        metadata = initialize_metadata(nb, pypackage=pypackage).dict()
        nb.metadata["nbproject"] = metadata
        write_notebook(nb, notebook_path)
        logger.success("attached notebook id to ipynb file")
    else:
        update_notebook_metadata(nb, notebook_path)
    return None


def save(notebook_path: str) -> Optional[str]:
    try:
        import nbstripout  # noqa
        from nbproject.dev import (
            MetaContainer,
            MetaStore,
            check_consecutiveness,
            read_notebook,
        )
        from nbproject.dev._meta_live import get_title
    except ImportError:
        logger.error("install nbproject & nbstripout: pip install nbproject nbstripout")
        return None
    nb = read_notebook(notebook_path)  # type: ignore
    nb_meta = nb.metadata
    is_consecutive = check_consecutiveness(nb)
    if not is_consecutive:
        if os.getenv("LAMIN_TESTING") is None:
            decide = input("   Do you still want to proceed with publishing? (y/n) ")
        else:
            decide = "n"
        if decide != "y":
            logger.error("Aborted (non-consecutive)!")
            return "aborted-non-consecutive"
    if get_title(nb) is None:
        logger.error(
            f"No title! Update & {colors.bold('save')} your notebook with a title '# My"
            " title' in the first cell."
        )
        return "no-title"
    if nb_meta is not None and "nbproject" in nb_meta:
        meta_container = MetaContainer(**nb_meta["nbproject"])
    else:
        logger.error("notebook isn't initialized, run lamin track <notebook_path>")
        return "not-initialized"

    meta_store = MetaStore(meta_container, notebook_path)
    import lamindb as ln

    ln.settings.verbosity = "success"
    transform_version = meta_store.version
    # the corresponding transform family in the transform table
    transform_family = ln.Transform.filter(uid__startswith=meta_store.id).all()
    if len(transform_family) == 0:
        logger.error(
            f"Didn't find notebook with stem_id {meta_store.id} (12 initial characters)"
            " in transform registry, did you run ln.track() in it?"
        )
        return "not-tracked-in-transform-registry"
    # the specific version
    transform = transform_family.filter(version=transform_version).one()
    # latest run of this transform by user
    run = ln.Run.filter(transform=transform).order_by("-run_at").first()
    if run.created_by.id != lamindb_setup.settings.user.id:
        response = input(
            "You are trying to save a notebook created by another user: Source and"
            " report files will be tagged with *your* user id. Proceed? (y/n)"
        )
        if response != "y":
            return "aborted-save-notebook-created-by-different-user"
    # convert the notebook file to html
    notebook_path_html = notebook_path.replace(".ipynb", ".html")
    logger.info(f"exporting notebook as html {notebook_path_html}")
    result = subprocess.run(f"jupyter nbconvert --to html {notebook_path}", shell=True)
    assert result.returncode == 0
    # copy the notebook file to a temporary file
    notebook_path_stripped = notebook_path.replace(".ipynb", "_stripped.ipynb")
    shutil.copy2(notebook_path, notebook_path_stripped)
    logger.info("stripping output of {notebook_path_tmp}")
    result = subprocess.run(f"nbstripout {notebook_path_stripped}", shell=True)
    assert result.returncode == 0
    # register the html report
    initial_report = None
    initial_source = None
    if len(transform_family) > 0:
        for transform in transform_family.order_by("-created_at"):
            # check for id to avoid query
            if transform.latest_report_id is not None:
                # any previous latest report of this transform is OK!
                initial_report = transform.latest_report
            if transform.source_file_id is not None:
                # any previous source file id is OK!
                initial_source = transform.source_file
    if hasattr(
        ln.settings, "silence_file_run_transform_warning"
    ):  # can remove this from next release on
        ln.settings.silence_file_run_transform_warning = True
    # register the source code
    if transform.source_file is not None:
        # this if condition is relevant if we already wrote the source code for
        # this transform version previously
        source_file = ln.File(
            notebook_path_stripped,
            description=f"Source of transform {transform.uid}",
            visibility=1,
        )
        # if the hash of the notebook file matches, we're trying to overwrite
        if source_file._state.adding:
            # this part of the if condition is currently not tested in an auto-mated way
            # will add this soon!
            if os.getenv("LAMIN_TESTING") is None:
                response = input(
                    "You're trying to save new notebook source code with the same"
                    " version; do you want to replace the existing source code? (y/n)"
                )
            else:
                response = "y"
            if response == "y":
                source_file = ln.File.filter(
                    version=transform_version, initial_version=initial_source
                ).one_or_none()
                if source_file is None:
                    # account for the very first version in the family
                    source_file = ln.File.filter(
                        version=transform_version, id=initial_source.id  # type: ignore
                    ).one_or_none()
                source_file.replace(notebook_path_stripped)
                source_file.save()
            else:
                logger.warning(
                    "Please create a new version of the notebook via `lamin track"
                    " <notebook_path>` and re-run the notebook"
                )
                return "rerun-the-notebook"
    else:
        source_file = ln.File(
            notebook_path_stripped,
            description=f"Source of transform {transform.uid}",
            version=transform_version,
            is_new_version_of=initial_source,
            visibility=1,
        )
        source_file.save()
    # save report file
    report_file = ln.File(
        notebook_path_html,
        description=f"Report of transform {transform.uid}",
        is_new_version_of=initial_report,
        visibility=1,
    )
    report_file.save()
    run.report = report_file
    run.save()
    # annotate transform
    transform.source_file = source_file
    transform.latest_report = report_file
    transform.save()
    # clean up
    Path(notebook_path_stripped).unlink()
    Path(notebook_path_html).unlink()
    msg = "saved notebook and wrote source file and html report"
    msg += (
        f"\n{transform}\ntransform.source_file: {source_file}\ntransform.latest_report:"
        f" {report_file}"
    )
    logger.success(msg)
    return None
