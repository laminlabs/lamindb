import os
import shutil
from pathlib import Path

import nox
from laminci import move_built_docs_to_docs_slash_project_slug, upload_docs_artifact
from lndb.test.nox import (
    build_docs,
    login_testuser1,
    login_testuser2,
    run_pre_commit,
    run_pytest,
)

import lamindb as ln

nox.options.reuse_existing_virtualenvs = True


@nox.session
def lint(session: nox.Session) -> None:
    run_pre_commit(session)


@nox.session
@nox.parametrize("package", ["bionty", "lnschema-bionty"])
def build(session, package):
    login_testuser2(session)
    login_testuser1(session)
    # run with pypi install on main
    if "GITHUB_EVENT_NAME" in os.environ:
        if os.environ["GITHUB_EVENT_NAME"] != "push":
            # run with submodule install on a PR
            session.install("./sub/lnschema-core[dev,test]")
            session.install("./sub/lnschema-wetlab[dev,test]")
            session.install("./sub/lndb-storage[dev,test]")
    session.install(".[dev,test]")

    session.install(".[dev,test]")
    if package == "lamindb":
        run_pytest(session)
    else:
        # navigate into submodule so that lamin-project.yml is correctly read
        os.chdir("./lnschema-bionty")
        session.install(".[test]")
        session.run("pytest", "-s", "./tests", "--ignore", "./tests/test_migrations.py")

    if package == "lamindb":
        # Schemas
        ln.setup.load("testuser1/lamin-site-assets")

        file = ln.select(ln.File, name="lnschema_core_docs").one()
        shutil.unpack_archive(file.load(), "lnschema_core_docs")
        Path("lnschema_core_docs/guide/0-core-schema.ipynb").rename(
            "docs/guide/lnschema-core.ipynb"
        )
        Path("lnschema_core_docs/guide/1-data-validation.ipynb").rename(
            "docs/guide/data-validation.ipynb"
        )

        file = ln.select(ln.File, name="lnschema_bionty_docs").one()
        shutil.unpack_archive(file.load(), "lnschema_bionty_docs")
        Path("lnschema_bionty_docs/guide/orms.ipynb").rename(
            "docs/guide/lnschema-bionty.ipynb"
        )
        Path("lnschema_bionty_docs/guide/knowledge.ipynb").rename(
            "docs/guide/knowledge.ipynb"
        )

        build_docs(session)
        login_testuser1(session)
        upload_docs_artifact()
        move_built_docs_to_docs_slash_project_slug()
