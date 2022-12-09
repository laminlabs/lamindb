from pathlib import Path

import nox
from lndb_setup.test.nox import build_docs, login_testuser1, run_pre_commit, run_pytest

nox.options.reuse_existing_virtualenvs = True


@nox.session(python=["3.7", "3.8", "3.9", "3.10", "3.11"])
def lint(session: nox.Session) -> None:
    run_pre_commit(session)


@nox.session(python=["3.7", "3.8", "3.9", "3.10", "3.11"])
def build(session):
    login_testuser1(session)
    login_user_2 = "lndb login testuser2@lamin.ai --password goeoNJKE61ygbz1vhaCVynGERaRrlviPBVQsjkhz"  # noqa
    session.run(*(login_user_2.split(" ")), external=True)
    session.install(".[dev,test]")
    instance_dirs = [
        d for d in ["./docs/guide/mydata", "./mydata-test-db"] if Path(d).exists()
    ]
    for instance_dir in instance_dirs:
        clean_instance = f"rm -r {instance_dir}"
        session.run(*clean_instance.split(" "))
    test_db = "lndb init --storage mydata-test-db"
    session.run(*test_db.split(" "), external=True)
    run_pytest(session)
    build_docs(session)
