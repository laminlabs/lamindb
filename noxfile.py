import os
import shutil
from pathlib import Path

import nox
from laminci import upload_docs_artifact
from laminci.nox import (
    build_docs,
    login_testuser1,
    login_testuser2,
    run,
    run_pre_commit,
)

# we'd like to aggregate coverage information across sessions
# and for this the code needs to be located in the same
# directory in every github action runner
# this also allows to break out an installation section
nox.options.default_venv_backend = "none"

IS_PR = os.getenv("GITHUB_EVENT_NAME") != "push"


GROUPS = {}
GROUPS["tutorial"] = [
    "introduction.ipynb",
    "tutorial.ipynb",
    "tutorial2.ipynb",
    "transfer.ipynb",
]
GROUPS["guide"] = [
    "data.ipynb",
    "meta.ipynb",
    "track.ipynb",
    "annotate.ipynb",
    "can-validate.ipynb",
    "annotate-for-developers.ipynb",
    "schemas.ipynb",
    "setup.ipynb",
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
    # on the release branch, do not use submodules but run with pypi install
    # only exception is the docs group
    if IS_PR or group == "docs":
        cmd = "uv pip install --system --no-deps ./sub/lamindb-setup ./sub/lnschema-core ./sub/lamin-cli"
        run(session, cmd)
    extras = ""
    if group == "unit":
        extras += "bionty,aws,zarr,fcs,jupyter"
    elif group == "tutorial":
        extras += "aws,jupyter,bionty"
    elif group == "guide":
        extras += "aws,bionty,zarr,jupyter,erdiagram"
        run(session, "uv pip install --system scanpy")
    elif group == "biology":
        extras += "bionty,fcs,jupyter"
    elif group == "faq":
        extras += "aws,bionty,jupyter"
    elif group == "storage":
        extras += "aws,zarr,bionty,jupyter"
        run(
            session,
            "uv pip install --system --no-deps ./sub/wetlab",
        )
        run(session, "uv pip install --system vitessce")
    elif group == "docs":
        extras += "bionty"
        run(session, "uv pip install --system mudata")
        run(
            session,
            "uv pip install --system --no-deps ./sub/wetlab",
        )
    elif group == "cli":
        extras += "jupyter,aws,bionty"
    if IS_PR and "bionty" in extras:
        run(
            session,
            "uv pip install --system --no-deps ./sub/bionty ./sub/lnschema-bionty",
        )
    run(session, f"uv pip install --system -e .[dev,{extras}]")


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
        run(session, f"pytest {coverage_args} ./tests")
    elif group == "tutorial":
        run(session, "lamin logout")
        run(
            session, f"pytest -s {coverage_args} ./docs/test_notebooks.py::test_{group}"
        )
    elif group == "guide":
        ln.setup.settings.auto_connect = True
        run(
            session,
            f"pytest -s {coverage_args} ./docs/test_notebooks.py::test_{group}",
        )
    elif group == "biology":
        run(
            session,
            f"pytest -s {coverage_args} ./docs/test_notebooks.py::test_{group}",
        )
    elif group == "faq":
        ln.setup.settings.auto_connect = True
        run(session, f"pytest -s {coverage_args} ./docs/faq")
    elif group == "storage":
        run(session, f"pytest -s {coverage_args} ./docs/storage")
    elif group == "cli":
        run(session, f"pytest {coverage_args} ./sub/lamin-cli/tests")
    # move artifacts into right place
    if group in {"tutorial", "guide", "biology"}:
        target_dir = Path(f"./docs/{group}")
        target_dir.mkdir(exist_ok=True)
        for filename in GROUPS[group]:
            shutil.copy(Path("docs") / filename, target_dir / filename)


@nox.session
def docs(session):
    # move artifacts into right place
    # for group in ["tutorial", "guide", "biology", "faq", "storage"]:
    #     if Path(f"./docs-{group}").exists():
    #         if Path(f"./docs/{group}").exists():
    #             shutil.rmtree(f"./docs/{group}")
    #         Path(f"./docs-{group}").rename(f"./docs/{group}")
    #     # move back to root level
    #     if group in {"tutorial", "guide", "biology"}:
    #         for path in Path(f"./docs/{group}").glob("*"):
    #             path.rename(f"./docs/{path.name}")
    run(session, "lamin init --storage ./docsbuild --schema bionty,wetlab")

    def generate_cli_docs():
        os.environ["NO_RICH"] = "1"
        from lamin_cli.__main__ import _generate_help

        page = "# `lamin`\n\nFor a guide, see: {doc}`/setup`.\n\n"
        helps = _generate_help()
        for name, help_string in helps.items():
            names = name.split(" ")
            section = ""
            if len(names) != 1:
                section = (
                    "```\n\n" + "#" * len(names) + " " + " ".join(("lamin", *names[1:]))
                )
            help_string = help_string.replace("Usage: main", "Usage: lamin")
            page += f"{section}\n\n```\n{help_string}```\n\n"
        Path("./docs/cli.md").write_text(page)

    generate_cli_docs()
    build_docs(session, strip_prefix=True, strict=True)
    upload_docs_artifact(aws=True)
