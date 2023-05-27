import os  # noqa
import shutil
from pathlib import Path
from subprocess import run
from time import perf_counter
from urllib.request import urlretrieve

import nox
from laminci import (  # noqa
    move_built_docs_to_docs_slash_project_slug,
    upload_docs_artifact,
)
from laminci.nox import login_testuser2  # noqa
from laminci.nox import build_docs, login_testuser1, run_pre_commit, run_pytest  # noqa

nox.options.reuse_existing_virtualenvs = True


@nox.session
def lint(session: nox.Session) -> None:
    run_pre_commit(session)


# we'd like to aggregate coverage information across sessions
# and for this the code needs to be located in the same
# directory in every github action runner
# hence, we're now running each session through a subprocess
# until we find a better solution
@nox.session
@nox.parametrize(
    "group",
    ["unit", "guide", "biology", "faq", "lndb-storage"],
)
def build(session, group):
    t_start = perf_counter()
    # run with pypi install on main
    if (
        "GITHUB_EVENT_NAME" in os.environ and os.environ["GITHUB_EVENT_NAME"] != "push"
    ):  # noqa
        # run with submodule install on a PR
        run("pip install --no-deps ./sub/lndb-setup", shell=True)
        run("pip install --no-deps ./sub/lnschema-core", shell=True)
        run("pip install --no-deps ./sub/lnbase-biolab", shell=True)
        run("pip install --no-deps ./sub/lndb-storage[dev,test]", shell=True)
    run("pip install .[dev,test]", shell=True)
    t_end = perf_counter()
    print(f"Done installing: {t_end - t_start:.3f}s")
    t_start = t_end

    login_testuser2(session)
    login_testuser1(session)
    t_end = perf_counter()
    print(f"Done logging in: {t_end - t_start:.3f}s")
    t_start = t_end

    if group != "lndb-storage":
        coverage_args = "--cov=lamindb --cov-append --cov-report=term-missing"  # noqa
        if group == "unit":
            run(f"pytest -s {coverage_args} ./tests", shell=True)
        elif group == "guide":
            run(f"pytest -s {coverage_args} ./docs/guide", shell=True)
        elif group == "biology":
            run(f"pytest -s {coverage_args} ./docs/biology", shell=True)
        elif group == "faq":
            run(f"pytest -s {coverage_args} ./docs/faq", shell=True)
    else:
        # navigate into submodule so that lamin-project.yml is correctly read
        os.chdir(f"./sub/{group}")
        run("pytest -s ./tests", shell=True)

    t_end = perf_counter()
    print(f"Done running tests: {t_end - t_start:.3f}s")


@nox.session
def docs(session):
    t_start = perf_counter()
    filename = "lndb_storage_docs.zip"
    urlretrieve(f"https://lamin-site-assets.s3.amazonaws.com/docs/{filename}", filename)
    shutil.unpack_archive(filename, "lndb_storage_docs")
    Path("lndb_storage_docs/guide/stream.ipynb").rename("docs/guide/stream.ipynb")

    filename = "lnschema_core_docs.zip"
    urlretrieve(f"https://lamin-site-assets.s3.amazonaws.com/docs/{filename}", filename)
    shutil.unpack_archive(filename, "lnschema_core_docs")
    Path("lnschema_core_docs/guide/0-core-schema.ipynb").rename(
        "docs/guide/lnschema-core.ipynb"
    )
    Path("lnschema_core_docs/guide/1-data-validation.ipynb").rename(
        "docs/guide/data-validation.ipynb"
    )

    filename = "lnschema_bionty_docs.zip"
    urlretrieve(f"https://lamin-site-assets.s3.amazonaws.com/docs/{filename}", filename)
    shutil.unpack_archive(filename, "lnschema_bionty_docs")
    Path("lnschema_bionty_docs/guide/bionty-orms.ipynb").rename(
        "docs/guide/lnschema-bionty.ipynb"
    )
    t_end = perf_counter()
    print(f"Done pulling artifacts: {t_end - t_start:.3f}s")
    t_start = t_end

    prefix = "." if Path("./lndocs").exists() else ".."
    run(f"pip install {prefix}/lndocs", shell=True)
    run("lamin init --storage ./docsbuild", shell=True)
    run("lndocs", shell=True)
    login_testuser1(session)
    upload_docs_artifact()
    move_built_docs_to_docs_slash_project_slug()

    t_end = perf_counter()
    print(f"Done building docs and uploading artifacts: {t_end - t_start:.3f}s")
