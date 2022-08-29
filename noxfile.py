from pathlib import Path

import nox

nox.options.reuse_existing_virtualenvs = True


@nox.session
def lint(session: nox.Session) -> None:
    session.install("pre-commit")
    session.run("pre-commit", "install")
    session.run("pre-commit", "run", "--all-files")


@nox.session(python=["3.9"])
def build(session):
    session.install(".[dev,test]")
    session.install(
        "git+https://github.com/tiangolo/sqlmodel.git"
    )  # just temporarily until the new release comes out to get rid of warnings # noqa
    login_user_1 = "lndb login raspbear@gmx.de --password MmR4YuQEyb0yxu7dAwJZTjLzR1Az2lN4Q4IduDlO"  # noqa
    login_user_2 = "lndb login kurt.hein@gmx.de --password uIoEGyiCj0qcXbGhTpOAuY6CH86xauzAsOSlp95A"  # noqa
    session.run(*(login_user_1.split(" ")))
    session.run(*(login_user_2.split(" ")))
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
