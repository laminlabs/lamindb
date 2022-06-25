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
        "lamindb",
        "setup",
        # "--notion",
        # os.environ["NOTION_API_KEY"],
        "--storage",
        "$HOME/data",
        "--user",
        "falexwolf",
    )
    session.run(
        "pytest",
        "--nbmake",
        "--overwrite",
    )  # write output instead of capturing it (more verbose)
    prefix = "." if Path("./lndocs").exists() else ".."
    session.install(f"{prefix}/lndocs")
    session.run("lndocs")
