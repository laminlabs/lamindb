import os  # noqa
import shutil
from pathlib import Path

import nox
from laminci import upload_docs_artifact
from laminci.nox import build_docs, login_testuser1, login_testuser2, run_pre_commit

# we'd like to aggregate coverage information across sessions
# and for this the code needs to be located in the same
# directory in every github action runner
# this also allows to break out an installation section
nox.options.default_venv_backend = "none"


GROUPS = {}
GROUPS["tutorial"] = [
    "introduction.ipynb",
    "tutorial.ipynb",
    "tutorial2.ipynb",
]
GROUPS["guide"] = [
    "data.ipynb",
    "meta.ipynb",
    "track.ipynb",
    "validate.ipynb",
    "annotate.ipynb",
    "schemas.ipynb",
    "setup.ipynb",
    "transfer.ipynb",
]
GROUPS["biology"] = [
    "public-ontologies.ipynb",
    "bio-registries.ipynb",
]


@nox.session
def lint(session: nox.Session) -> None:
    run_pre_commit(session)


@nox.session
@nox.parametrize(
    "group",
    ["unit", "tutorial", "guide", "biology", "faq", "storage", "docs", "cli"],
)
def install(session, group):
    # run with pypi install on main
    if os.getenv("GITHUB_EVENT_NAME") != "push":
        # run with submodule install on a PR
        submodules = " ".join(
            [
                "./sub/lamindb-setup",
                "./sub/lnschema-core",
                "./sub/lamin-cli",
            ]
        )
        session.run(*f"pip install --no-deps {submodules}".split())
    extras = ""
    if group == "unit":
        extras += "bionty,aws,zarr,postgres,fcs"
    elif group == "tutorial":
        extras += "aws,jupyter,bionty"  # despite no AWS credentials, we need s3fs
    elif group == "guide":
        extras += "aws,bionty,zarr,jupyter,erdiagram,postgres"
        session.run(*"pip install scanpy".split())
    elif group == "biology":
        extras += "bionty,fcs,jupyter"
    elif group == "faq":
        extras += "aws,postgres,bionty,jupyter"
    elif group == "storage":
        extras += "aws,zarr"
    elif group == "docs":
        extras += "bionty"
    elif group == "cli":
        extras += "jupyter"
    if os.getenv("GITHUB_EVENT_NAME") != "push":
        if "bionty" in extras:
            session.run(*"pip install --no-deps ./sub/lnschema-bionty".split())
    session.run(*f"pip install -e .[dev,{extras}]".split())


@nox.session
@nox.parametrize(
    "group",
    ["unit", "tutorial", "guide", "biology", "faq", "storage", "cli"],
)
def build(session, group):
    login_testuser2(session)
    login_testuser1(session)
    coverage_args = "--cov=lamindb --cov-append --cov-report=term-missing"
    if group == "unit":
        session.run(*f"pytest {coverage_args} ./tests".split())
    elif group == "tutorial":
        session.run(
            *f"pytest -s {coverage_args} ./docs/test_notebooks.py::test_{group}".split()
        )
    elif group == "guide":
        session.run(
            *f"pytest -s {coverage_args} ./docs/test_notebooks.py::test_{group}".split()
        )
    elif group == "biology":
        session.run(
            *f"pytest -s {coverage_args} ./docs/test_notebooks.py::test_{group}".split()
        )
    elif group == "faq":
        session.run(*f"pytest -s {coverage_args} ./docs/faq".split())
    elif group == "storage":
        session.run(*f"pytest -s {coverage_args} ./docs/storage".split())
    elif group == "cli":
        session.run(*f"pytest {coverage_args} ./sub/lamin-cli/tests".split())
    # move artifacts into right place
    if group in {"tutorial", "guide", "biology"}:
        target_dir = Path(f"./docs/{group}")
        target_dir.mkdir(exist_ok=True)
        for filename in GROUPS[group]:
            shutil.copy(Path("docs") / filename, target_dir / filename)


@nox.session
def docs(session):
    # move artifacts into right place
    for group in ["tutorial", "guide", "biology", "faq", "storage"]:
        if Path(f"./docs-{group}").exists():
            if Path(f"./docs/{group}").exists():
                shutil.rmtree(f"./docs/{group}")
            Path(f"./docs-{group}").rename(f"./docs/{group}")
        # move back to root level
        if group in {"tutorial", "guide", "biology"}:
            for path in Path(f"./docs/{group}").glob("*"):
                path.rename(f"./docs/{path.name}")
    login_testuser1(session)
    session.run(*"lamin init --storage ./docsbuild --schema bionty".split())

    def generate_cli_docs(main_parser):
        page = "# `lamin`\n\nFor a guide, see: {doc}`/setup`.\n\n"
        commands = [
            "login",
            "init",
            "load",
            "close",
            "delete",
            "track",
            "info",
            "migrate",
            "save",
            "set",
            "schema",
        ]
        for action_group in main_parser._action_groups:
            for group_action in action_group._group_actions:
                if type(group_action).__name__ == "_SubParsersAction":
                    for command in commands:
                        subparser = group_action.choices[command]
                        # replace the "nox" command with the "lamin" command
                        help_string = subparser.format_help().replace("nox", "lamin")
                        page += f"## `lamin {command}`\n\n```\n{help_string}```\n\n"
        Path("./docs/cli.md").write_text(page)

    from lamin_cli import __main__

    generate_cli_docs(__main__.parser)

    build_docs(session, strip_prefix=True, strict=True)
    upload_docs_artifact(aws=True)
