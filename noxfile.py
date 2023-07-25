import os  # noqa
import shutil
from pathlib import Path

import nox
from laminci import upload_docs_artifact
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
    ["unit", "guide", "biology", "faq", "storage", "docs"],
)
def install(session, group):
    # run with pypi install on main
    if os.getenv("GITHUB_EVENT_NAME") != "push":
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
        extras += ",bionty,aws,zarr,postgres"
    elif group == "guide":
        extras += ",aws,bionty,zarr,jupyter"
        session.run(*"pip install scanpy".split())
    elif group == "biology":
        extras += ",bionty,fcs,jupyter"
        session.run(*"pip install scanpy".split())
        session.run(*"pip install mudata".split())
    elif group == "faq":
        extras += ",aws,postgres,bionty,jupyter"
    elif group == "storage":
        extras += ",aws,zarr"
    elif group == "docs":
        extras += ",bionty"
    if os.getenv("GITHUB_EVENT_NAME") != "push":
        if "bionty" in extras:
            session.run(*"pip install --no-deps ./sub/lnschema-bionty".split())
    session.run(*f"pip install -e .[dev{extras}]".split())


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
        session.run(*f"pytest {coverage_args} ./tests".split())
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
        if group == "biology":
            Path("./docs/biology/lnschema-bionty.ipynb").rename(
                "./docs/lnschema-bionty.ipynb"
            )
    login_testuser1(session)
    session.run(*"lamin init --storage ./docsbuild --schema bionty".split())
    build_docs(session, strip_prefix=True, strict=True)
    upload_docs_artifact(aws=True)
