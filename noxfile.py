import os

import nox
from laminci import upload_docs_dir
from lndb.test.nox import (
    build_docs,
    login_testuser1,
    login_testuser2,
    run_pre_commit,
    run_pytest,
)

nox.options.reuse_existing_virtualenvs = True


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
    build_docs(session)
    login_testuser1(session)
    upload_docs_dir()
