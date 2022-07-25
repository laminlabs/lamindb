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
    # the following is necessary to make the CLI available in the virtual env
    # for the user, a typical pip install lamindb will also make the lndb
    # CLI available without the dedicated pip install of it!
    # Actually, probably it wasn't necessary and just clearing the cache was
    session.install("lndb_cli")
    session.install(".[dev,test]")
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
