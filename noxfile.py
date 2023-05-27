import os  # noqa
import shutil
from pathlib import Path
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


@nox.session
@nox.parametrize("package", ["lamindb", "lndb-storage"])
def build(session, package):
    t_start = perf_counter()
    # run with pypi install on main
    if (
        "GITHUB_EVENT_NAME" in os.environ and os.environ["GITHUB_EVENT_NAME"] != "push"
    ):  # noqa
        # run with submodule install on a PR
        session.install("./sub/lndb-setup")
        session.install("./sub/lnschema-core")
        session.install("./sub/lnbase-biolab")
        session.install("./sub/lndb-storage[dev,test]")
    session.install(".[dev,test]")
    t_end = perf_counter()
    print(f"Done installing: {t_end - t_start:.3f}s")
    t_start = t_end

    login_testuser2(session)
    login_testuser1(session)
    t_end = perf_counter()
    print(f"Done logging in: {t_end - t_start:.3f}s")
    t_start = t_end

    if package == "lamindb":
        run_pytest(session)
    else:
        # navigate into submodule so that lamin-project.yml is correctly read
        os.chdir(f"./sub/{package}")
        session.run(
            "pytest", "-s", "./tests", "--ignore", "./tests/test_migrations.py"
        )  # noqa

    t_end = perf_counter()
    # print(f"Done running tests: {t_end - t_start:.3f}s")
    t_start = t_end

    if package == "lamindb":
        filename = "lndb_storage_docs.zip"
        urlretrieve(
            f"https://lamin-site-assets.s3.amazonaws.com/docs/{filename}", filename
        )
        shutil.unpack_archive(filename, "lndb_storage_docs")
        Path("lndb_storage_docs/guide/stream.ipynb").rename("docs/guide/stream.ipynb")

        filename = "lnschema_core_docs.zip"
        urlretrieve(
            f"https://lamin-site-assets.s3.amazonaws.com/docs/{filename}", filename
        )
        shutil.unpack_archive(filename, "lnschema_core_docs")
        Path("lnschema_core_docs/guide/0-core-schema.ipynb").rename(
            "docs/guide/lnschema-core.ipynb"
        )
        Path("lnschema_core_docs/guide/1-data-validation.ipynb").rename(
            "docs/guide/data-validation.ipynb"
        )

        filename = "lnschema_bionty_docs.zip"
        urlretrieve(
            f"https://lamin-site-assets.s3.amazonaws.com/docs/{filename}", filename
        )
        shutil.unpack_archive(filename, "lnschema_bionty_docs")
        Path("lnschema_bionty_docs/guide/bionty-orms.ipynb").rename(
            "docs/guide/lnschema-bionty.ipynb"
        )
        t_end = perf_counter()
        print(f"Done pulling artifacts: {t_end - t_start:.3f}s")
        t_start = t_end

        build_docs(session)
        login_testuser1(session)
        upload_docs_artifact()
        move_built_docs_to_docs_slash_project_slug()

        t_end = perf_counter()
        print(f"Done building docs and uploading artifacts: {t_end - t_start:.3f}s")
