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
CI = os.environ.get("CI")


GROUPS = {}
GROUPS["tutorial"] = [
    "transfer.ipynb",
    "arrays.ipynb",
    "registries.ipynb",
    "setup.ipynb",
]
GROUPS["guide"] = [
    "track.ipynb",
    "curate-df.ipynb",
    "curate-any.ipynb",
    "curate-subclass.ipynb",
    "schemas.ipynb",
]
GROUPS["biology"] = [
    "bio-registries.ipynb",
]


@nox.session
def lint(session: nox.Session) -> None:
    run_pre_commit(session)


@nox.session
def install(session):
    base_deps = [
        "./sub/lnschema-core",
        "./sub/lamin-cli",
        "./sub/lamindb-setup",
        "./sub/bionty",
    ]
    top_deps = [
        ".[aws,bionty,jupyter]",
    ]
    cmds = [
        f"uv pip install {'--system' if CI else ''} --no-cache-dir {' '.join(base_deps)}",
    ] + [
        f"uv pip install {'--system' if CI else ''} --no-cache-dir -e {dep}"
        for dep in top_deps
    ]
    [run(session, line) for line in cmds]


@nox.session
@nox.parametrize(
    "group",
    [
        "unit-core",
        "unit-storage",
        "tutorial",
        "guide",
        "biology",
        "faq",
        "storage",
        "spatial",
        "docs",
        "cli",
    ],
)
def install_ci(session, group):
    extras = ""
    if group == "unit-core":
        extras += "bionty,aws,gcp,zarr,fcs,jupyter"
        # testing load_to_memory for yaml
        run(session, "uv pip install --system PyYAML")
        run(session, "uv pip install --system huggingface_hub")
        # pinning 1.15.0rc3 because 1.14.5 is incompatible with anndata>=0.11.0
        # and >1.15.0rc4 has a different append mode API
        # wating for a normal release to adapt
        run(
            session, "uv pip install --system tiledbsoma==1.15.0rc3"
        )  # test SOMACurator
    elif group == "unit-storage":
        extras += "aws,zarr,bionty"
        run(session, "uv pip install --system tiledbsoma==1.15.0rc3")
    elif group == "tutorial":
        extras += "aws,jupyter,bionty"
        run(session, "uv pip install --system huggingface_hub")
    elif group == "guide":
        extras += "aws,bionty,zarr,jupyter"
        run(session, "uv pip install --system scanpy")
    elif group == "biology":
        extras += "bionty,fcs,jupyter"
        run(session, "uv pip install --system ipywidgets")
    elif group == "faq":
        extras += "aws,bionty,jupyter"
        run(
            session,
            "uv pip install --system --no-deps ./sub/findrefs ./sub/ourprojects",
        )
    elif group == "storage":
        extras += "aws,zarr,bionty,jupyter"
        run(
            session,
            "uv pip install --system --no-deps ./sub/wetlab ./sub/findrefs ./sub/ourprojects",
        )
        run(session, "uv pip install --system vitessce")
    elif group == "spatial":
        extras += "aws,zarr,bionty,jupyter"
        run(
            session,
            "uv pip install --system ./sub/lamin-spatial",
        )
        run(
            session,
            "uv pip install --system -U git+https://github.com/scverse/spatialdata.git@refs/pull/806/head",
        )  # Required to access metadata attrs
    elif group == "docs":
        extras += "bionty,zarr"
        run(session, "uv pip install --system mudata")
        run(
            session,
            "uv pip install --system --no-deps ./sub/wetlab ./sub/findrefs ./sub/clinicore ./sub/omop ./sub/cellregistry ./sub/ourprojects .sub/lamin-spatial",
        )
    elif group == "cli":
        extras += "jupyter,aws,bionty"
    run(session, f"uv pip install --system -e .[dev,{extras}]")
    # on the release branch, do not use submodules but run with pypi install
    # only exception is the docs group which should always use the submodule
    # to push docs fixes fast
    # installing this after lamindb to be sure that these packages won't be reinstaled
    # during lamindb installation
    if IS_PR or group == "docs":
        cmd = "uv pip install --system --no-deps ./sub/lamindb-setup ./sub/lnschema-core ./sub/lamin-cli"
        run(session, cmd)
        if "bionty" in extras:
            run(
                session,
                "uv pip install --system --no-deps ./sub/bionty",
            )


@nox.session
@nox.parametrize(
    "group",
    [
        "unit-core",
        "unit-storage",
        "spatial",
        "tutorial",
        "guide",
        "biology",
        "faq",
        "storage",
        "cli",
    ],
)
def build(session, group):
    import lamindb as ln

    login_testuser2(session)
    login_testuser1(session)
    run(session, "lamin settings set private-django-api true")
    coverage_args = "--cov=lamindb --cov-config=pyproject.toml --cov-append --cov-report=term-missing"
    if group == "unit-core":
        run(session, f"pytest {coverage_args} ./tests/core --durations=50")
    elif group == "unit-storage":
        run(session, f"pytest {coverage_args} ./tests/storage --durations=50")
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
    elif group == "spatial":
        run(
            session,
            f"pytest {coverage_args} ./sub/lamin-spatial/tests/test_spatialdata_curator.py --durations=50",
        )
    elif group == "cli":
        run(session, f"pytest {coverage_args} ./sub/lamin-cli/tests --durations=50")
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
    run(
        session,
        "lamin init --storage ./docsbuild --schema bionty,wetlab,clinicore,findrefs,ourprojects,cellregistry,omop",
    )

    def generate_cli_docs():
        os.environ["NO_RICH"] = "1"
        from lamin_cli.__main__ import _generate_help

        page = "# `lamin`\n\nFor a guide, see: {doc}`/setup`.\n"
        helps = _generate_help()
        for name, help_dict in helps.items():
            names = name.split(" ")
            section = ""
            if len(names) != 1:
                if len(names) == 2:
                    separator = "#" * len(names) + " " + " ".join(("lamin", *names[1:]))
                else:
                    separator = ""
                section = "```\n\n" + separator

            help_string = help_dict["help"].replace("Usage: main", "Usage: lamin")
            help_docstring = help_dict["docstring"]
            page += f"{section}\n\n{help_docstring}\n\n```text\n{help_string}\n\n"
        Path("./docs/cli.md").write_text(page)

    generate_cli_docs()
    build_docs(session, strip_prefix=True, strict=True)
    upload_docs_artifact(aws=True)
