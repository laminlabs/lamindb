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
    t_end = perf_counter()
    print(f"Done logging in: {t_end - t_start:.3f}s")
    t_start = t_end

    # run with pypi install on main
    if "GITHUB_EVENT_NAME" in os.environ and os.environ["GITHUB_EVENT_NAME"] != "push":
        # run with submodule install on a PR
        session.install("./sub/lndb-setup")
        session.install("./sub/lnschema-core")
        session.install("./sub/lnbase-biolab")
        session.install("./sub/lndb-storage[dev,test]")
    session.install(".[dev,test]")
    t_end = perf_counter()
    print(f"Done installing: {t_end - t_start:.3f}s")
    t_start = t_end

    if package == "lamindb":
        run_pytest(session)
    else:
        # navigate into submodule so that lamin-project.yml is correctly read
        os.chdir(f"./sub/{package}")
        session.run("pytest", "-s", "./tests", "--ignore", "./tests/test_migrations.py")

    t_end = perf_counter()
    print(f"Done running tests: {t_end - t_start:.3f}s")
    t_start = t_end

    if package == "lamindb":
        # Schemas
        ln.setup.load("testuser1/lamin-site-assets", migrate=True)

        file = ln.select(ln.File, key="docs/lndb_storage_docs.zip").one()
        shutil.unpack_archive(file.stage(), "lndb_storage_docs")
        Path("lndb_storage_docs/guide/stream.ipynb").rename("docs/guide/stream.ipynb")

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

        t_end = perf_counter()
        print(f"Done pulling artifacts: {t_end - t_start:.3f}s")
        t_start = t_end

        ln.setup.init(storage="./mydocs")  # instance for building docs
        build_docs(session)
        login_testuser1(session)
        upload_docs_artifact()
        move_built_docs_to_docs_slash_project_slug()

        t_end = perf_counter()
        print(f"Done building docs and uploading artifacts: {t_end - t_start:.3f}s")
