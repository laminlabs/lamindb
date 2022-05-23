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
    session.run(
        "pytest",
        "--nbmake",
        "--overwrite",
    )  # write output instead of capturing it (more verbose)
    if Path("./lndocs").exists():  # GitHub Actions prefers nesting
        session.install("./lndocs")
    else:
        session.install("../lndocs")  # Locally we have a flat structure
    session.install(
        "git+https://github.com/pydata/pydata-sphinx-theme.git"
    )  # just temporarily until the new release comes out
    session.run("lndocs")
