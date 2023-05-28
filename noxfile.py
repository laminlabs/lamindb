import os  # noqa
import shutil
from pathlib import Path
from urllib.request import urlretrieve

import nox
from laminci import (  # noqa
    move_built_docs_to_docs_slash_project_slug,
    upload_docs_artifact,
)
from laminci.nox import login_testuser2  # noqa
from laminci.nox import build_docs, login_testuser1, run_pre_commit, run_pytest  # noqa


@nox.session
def lint(session: nox.Session) -> None:
    session.run(*"pip install pre-commit".split())
    session.run("pre-commit", "install")
    session.run("pre-commit", "run", "--all-files")


# we'd like to aggregate coverage information across sessions
# and for this the code needs to be located in the same
# directory in every github action runner
# this also allows to break out an installation section
nox.options.default_venv_backend = "none"


@nox.session
def install(session):
    # run with pypi install on main
    if "GITHUB_EVENT_NAME" in os.environ and os.environ["GITHUB_EVENT_NAME"] != "push":
        # run with submodule install on a PR
        session.run(*"pip install --no-deps ./sub/lndb-setup".split())
        session.run(*"pip install --no-deps ./sub/lnschema-core".split())
        session.run(*"pip install --no-deps ./sub/lnbase-biolab".split())
        session.run(*"pip install --no-deps ./sub/lndb-storage".split())
    session.run(*"pip install .[dev,test]".split())


@nox.session
@nox.parametrize(
    "group",
    ["unit", "guide", "biology", "faq", "lndb-storage"],
)
def build(session, group):
    login_testuser2(session)
    login_testuser1(session)
    coverage_args = "--cov=lamindb --cov-append --cov-report=term-missing"  # noqa
    if group == "unit":
        session.run(*f"pytest -s {coverage_args} ./tests".split())
    elif group == "guide":
        session.run(*f"pytest -s {coverage_args} ./docs/guide".split())
    elif group == "biology":
        session.run(*f"pytest -s {coverage_args} ./docs/biology".split())
    elif group == "faq":
        session.run(*f"pytest -s {coverage_args} ./docs/faq".split())
    elif group == "lndb-storage":
        with session.chdir(f"./sub/{group}"):
            session.run(*f"pytest -s {coverage_args} ./tests".split())
            # I'd like to replace below with
            # Path(f"./sub/{group}/.coverage").rename(".")
            # but it errored with
            # OSError: Device or resource busy: 'sub/lndb-storage/.coverage' -> '.'
            session.run(*"mv .coverage ../..".split())


@nox.session
def docs(session):
    # move artifacts into right place
    for group in ["guide", "biology", "faq"]:
        if Path(f"./docs-{group}").exists():
            shutil.rmtree(f"./docs/{group}")
            Path(f"./docs-{group}").rename(f"./docs/{group}")

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

    prefix = "." if Path("./lndocs").exists() else ".."
    session.run(*f"pip install {prefix}/lndocs".split())
    login_testuser1(session)
    session.run(*"lamin init --storage ./docsbuild".split())
    session.run("lndocs")
    upload_docs_artifact()
    move_built_docs_to_docs_slash_project_slug()
