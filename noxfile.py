import os
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
        extras += "bionty,aws,zarr,postgres,fcs,jupyter"
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
        extras += "aws,zarr,bionty,jupyter,postgres"
        session.run(
            *"pip install --no-deps wetlab@git+https://github.com/laminlabs/wetlab".split()
        )
    elif group == "docs":
        extras += "bionty"
    elif group == "cli":
        extras += "jupyter,aws"
    if os.getenv("GITHUB_EVENT_NAME") != "push" and "bionty" in extras:
        session.run(*"pip install --no-deps ./sub/bionty".split())
        session.run(*"pip install --no-deps ./sub/lnschema-bionty".split())
    session.run(*f"pip install -e .[dev,{extras}]".split())


@nox.session
@nox.parametrize(
    "group",
    ["unit", "tutorial", "guide", "biology", "faq", "storage", "cli"],
)
def build(session, group):
    import lamindb as ln

    login_testuser2(session)
    login_testuser1(session)
    coverage_args = "--cov=lamindb --cov-append --cov-report=term-missing"
    if group == "unit":
        session.run(*f"pytest {coverage_args} ./tests".split())
    elif group == "tutorial":
        session.run(*"lamin logout".split())
        session.run(
            *f"pytest -s {coverage_args} ./docs/test_notebooks.py::test_{group}".split()
        )
    elif group == "guide":
        ln.setup.settings.auto_connect = True
        session.run(
            *f"pytest -s {coverage_args} ./docs/test_notebooks.py::test_{group}".split()
        )
    elif group == "biology":
        session.run(
            *f"pytest -s {coverage_args} ./docs/test_notebooks.py::test_{group}".split()
        )
    elif group == "faq":
        ln.setup.settings.auto_connect = True
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
    session.run(*"lamin init --storage ./docsbuild --schema bionty".split())

    def generate_cli_docs():
        os.environ["NO_RICH"] = "1"
        from lamin_cli.__main__ import _generate_help

        page = "# `lamin`\n\nFor a guide, see: {doc}`/setup`.\n\n"
        helps = _generate_help()

        for name, help_string in helps.items():
            print(name, help_string)
            names = name.split(" ")
            section = ""
            if len(names) != 1:
                section = (
                    "```\n\n" + "#" * len(names) + " " + " ".join(("lamin", *names[1:]))
                )
            page += f"{section}\n\n```\n{help_string}```\n\n"

        Path("./docs/cli.md").write_text(page)

    generate_cli_docs()
    build_docs(session, strip_prefix=True, strict=True)
    upload_docs_artifact(aws=True)
