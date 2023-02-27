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
    session.install(".[dev,test]")
    run_pytest(session)
    build_docs(session)
    upload_docs_dir()
