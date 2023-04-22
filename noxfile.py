import os
import shutil
from pathlib import Path
from time import perf_counter

import nox
from laminci import move_built_docs_to_docs_slash_project_slug, upload_docs_artifact
from laminci.nox import (
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
@nox.parametrize("package", ["lamindb", "lndb-storage"])
def build(session, package):
    t_start = perf_counter()
    login_testuser2(session)
    login_testuser1(session)
    t_total = perf_counter() - t_start
    print(f"Done logging in: {t_total:.3f}s")

    t_start = perf_counter()
    # run with pypi install on main
    if "GITHUB_EVENT_NAME" in os.environ and os.environ["GITHUB_EVENT_NAME"] != "push":
        # run with submodule install on a PR
        session.install("./sub/lnschema-core[dev,test]")
        session.install("./sub/lnbase-biolab[dev,test]")
        session.install("./sub/lndb-storage[dev,test]")

    session.install(".[dev,test]")
    t_total = perf_counter() - t_start
    print(f"Done installing: {t_total:.3f}s")

    if package == "lamindb":
        run_pytest(session)
    else:
        # navigate into submodule so that lamin-project.yml is correctly read
        os.chdir(f"./sub/{package}")
        session.run("pytest", "-s", "./tests", "--ignore", "./tests/test_migrations.py")

    if package == "lamindb":
        t_start = perf_counter()

        # Schemas
        ln.setup.load("testuser1/lamin-site-assets", migrate=True)

        file = ln.select(ln.File, key="docs/lnschema_core_docs.zip").one()
        shutil.unpack_archive(file.stage(), "lnschema_core_docs")
        Path("lnschema_core_docs/guide/0-core-schema.ipynb").rename(
            "docs/guide/lnschema-core.ipynb"
        )
        Path("lnschema_core_docs/guide/1-data-validation.ipynb").rename(
            "docs/guide/data-validation.ipynb"
        )

        file = ln.select(ln.File, key="docs/lnschema_bionty_docs.zip").one()
        shutil.unpack_archive(file.stage(), "lnschema_bionty_docs")
        Path("lnschema_bionty_docs/guide/bionty-orms.ipynb").rename(
            "docs/guide/lnschema-bionty.ipynb"
        )
        Path("lnschema_bionty_docs/guide/knowledge.ipynb").rename(
            "docs/guide/knowledge.ipynb"
        )

        t_total = perf_counter() - t_start
        print(f"Done pulling artifacts: {t_total:.3f}s")

        t_start = perf_counter()
        build_docs(session)
        login_testuser1(session)
        upload_docs_artifact()
        move_built_docs_to_docs_slash_project_slug()
        print(f"Done building docs and uploading: {t_total:.3f}s")
