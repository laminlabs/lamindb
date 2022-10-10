from pathlib import Path

import nox

nox.options.reuse_existing_virtualenvs = True
nox.options.error_on_external_run = False
nox.options.default_venv_backend = None


@nox.session(python=["3.7", "3.8", "3.9", "3.10", "3.11"])
def lint(session: nox.Session) -> None:
    session.install("pre-commit")
    session.run("pre-commit", "install")
    session.run("pre-commit", "run", "--all-files")


@nox.session(python=["3.7", "3.8", "3.9", "3.10", "3.11"])
def build(session):
    session.install(".[dev,test]")
    login_user_1 = "lndb login testuser1@lamin.ai --password cEvcwMJFX4OwbsYVaMt2Os6GxxGgDUlBGILs2RyS"  # noqa
    login_user_2 = "lndb login testuser2@lamin.ai --password goeoNJKE61ygbz1vhaCVynGERaRrlviPBVQsjkhz"  # noqa
    session.run(*(login_user_1.split(" ")))
    session.run(*(login_user_2.split(" ")))
    test_db = "lndb init --storage mydata-test-db"
    session.run(*test_db.split(" "))
    session.run(
        "pytest",
        "-s",
        "--cov=lamindb",
        "--cov-append",
        "--cov-report=term-missing",
    )
    session.run("coverage", "xml")
    prefix = "." if Path("./lndocs").exists() else ".."
    session.install(f"{prefix}/lndocs")
    session.run("lndocs")
