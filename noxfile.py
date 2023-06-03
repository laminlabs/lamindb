import os  # noqa
import shutil
from pathlib import Path
from urllib.request import urlretrieve

import nox
from laminci import (  # noqa
    move_built_docs_to_docs_slash_project_slug,
    upload_docs_artifact,
)
from laminci.nox import build_docs, login_testuser1, login_testuser2, run_pre_commit

# we'd like to aggregate coverage information across sessions
# and for this the code needs to be located in the same
# directory in every github action runner
# this also allows to break out an installation section
nox.options.default_venv_backend = "none"


@nox.session
def lint(session: nox.Session) -> None:
    run_pre_commit(session)


@nox.session
@nox.parametrize(
    "group",
    ["unit", "guide", "biology", "faq", "storage"],
)
def install(session, group):
    # run with pypi install on main (currently disabled)
    if True:
        # os.getenv("GITHUB_EVENT_NAME") not in {None, "push"}
        # run with submodule install on a PR
        submodules = " ".join(
            [
                "./sub/lamindb-setup",
                "./sub/lnschema-core",
            ]
        )
        session.run(*f"pip install --no-deps {submodules}".split())
    extras = ""
    if group == "unit":
        extras += ",bionty"
    elif group == "guide":
        extras += ",aws"
    elif group == "biology":
        extras += ",lamin1"
    elif group == "storage":
        extras += ",aws"
    session.run(*f"pip install .[test{extras}]".split())


@nox.session
@nox.parametrize(
    "group",
    ["unit", "guide", "biology", "faq", "storage"],
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
    elif group == "storage":
        session.run(*f"pytest -s {coverage_args} ./docs/storage".split())


@nox.session
def docs(session):
    # move artifacts into right place
    for group in ["guide", "biology", "faq", "storage"]:
        if Path(f"./docs-{group}").exists():
            shutil.rmtree(f"./docs/{group}")
            Path(f"./docs-{group}").rename(f"./docs/{group}")

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

    login_testuser1(session)
    session.run(*"lamin init --storage ./docsbuild".split())
    session.run(*"pip install django".split())
    build_docs(session)
    # do not upload docs artifact for now
    # upload_docs_artifact()
    move_built_docs_to_docs_slash_project_slug()
