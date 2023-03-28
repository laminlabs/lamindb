import os
import shutil
from pathlib import Path
from typing import List, Tuple

import nox
from laminci import upload_docs_dir
from lndb.test.nox import (
    build_docs,
    login_testuser1,
    login_testuser2,
    run_pre_commit,
    run_pytest,
)

import lamindb as ln

nox.options.reuse_existing_virtualenvs = True


def replace_content(filename: Path, mapped_content: List[Tuple[str, str]]) -> None:
    with open(filename) as f:
        content = f.read()
    with open(filename, "w") as f:
        for args in mapped_content:
            content = content.replace(*args)
        f.write(content)


@nox.session(python=["3.7", "3.8", "3.9", "3.10", "3.11"])
def lint(session: nox.Session) -> None:
    run_pre_commit(session)


@nox.session(python=["3.7", "3.8", "3.9", "3.10", "3.11"])
def build(session):
    login_testuser1(session)
    login_testuser2(session)
    # run with pypi install on main
    if "GITHUB_EVENT_NAME" in os.environ:
        if os.environ["GITHUB_EVENT_NAME"] != "push":
            # run with submodule install on a PR
            session.install("./lnschema-core[dev,test]")
            session.install("./lnschema-wetlab[dev,test]")
            session.install("./lndb-storage[dev,test]")
    session.install(".[dev,test]")
    run_pytest(session)

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

    # Bionty

    file = ln.select(ln.File, name="bionty_docs").one()
    shutil.unpack_archive(file.load(), "bionty_docs")
    Path("bionty_docs").rename("docs/bionty")

    replace_content(
        "docs/bionty/index.md",
        mapped_content=[
            ("../README.md", "./README.md"),
        ],
    )

    build_docs(session)
    login_testuser1(session)
    upload_docs_dir()
