import nox

nox.options.reuse_existing_virtualenvs = True


@nox.session
def lint(session: nox.Session) -> None:
    session.install("pre-commit")
    session.run("pre-commit", "install")
    session.run("pre-commit", "run", "--all-files")


@nox.session
def build(session):
    session.install(".[dev, test]")
    session.run(
        "pytest",
        "--nbmake",
        "--overwrite",
    )  # write output instead of capturing it (more verbose)
    session.install("-r", "docs/lamin_sphinx/requirements.txt")
    session.run(*"sphinx-build -b html docs _build/html".split())
