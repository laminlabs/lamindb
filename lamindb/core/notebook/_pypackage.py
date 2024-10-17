from __future__ import annotations

import re
import sys
from ast import Import, ImportFrom, parse, walk
from typing import TYPE_CHECKING

from importlib_metadata import PackageNotFoundError, packages_distributions, version

if TYPE_CHECKING:
    from collections.abc import Iterable

    from nbformat import NotebookNode

std_libs = None
pkgs_dists = None


def _load_pkgs_info():
    global std_libs
    global pkgs_dists

    major, minor = sys.version_info[0], sys.version_info[1]
    if major == 3 and minor > 9:
        std_libs = sys.stdlib_module_names  # type: ignore
    else:
        from stdlib_list import stdlib_list

        std_libs = set(stdlib_list(f"{major}.{minor}"))

    pkgs_dists = packages_distributions()


def _get_version(pkg):
    try:
        pkg_ver = version(pkg)
    except PackageNotFoundError:
        pkg_ver = ""
    return pkg_ver


def cell_imports(cell_source: str):
    # based on the package https://github.com/bndr/pipreqs for python scripts
    # parses python import statements in the code cells
    tree = parse(cell_source)
    for node in walk(tree):
        if isinstance(node, Import):
            for subnode in node.names:
                name = subnode.name.partition(".")[0]
                if name != "":
                    yield name
        elif isinstance(node, ImportFrom):
            name = node.module.partition(".")[0]  # type: ignore
            if name != "":
                yield name


def infer_pypackages(
    content: NotebookNode,
    add_pkgs: Iterable | None = None,
    pin_versions: bool = True,
):
    """Parse notebook object and infer all pypackages.

    This does not account for additional packages in file metadata.

    For the user-facing functionality,
    see :meth:`~nbproject.dev.MetaLive.pypackage`.

    Args:
        content: A notebook to infer pypackages from.
        add_pkgs: Additional packages to add.
        pin_versions: If `True`, fixes versions from the current environment.

    Examples:
        >>> pypackages = nbproject.dev.infer_pypackages(nb)
        >>> pypackages
        {"scanpy": "1.8.7", "pandas": "1.4.3"}
    """
    cells = content.cells

    if std_libs is None or pkgs_dists is None:
        _load_pkgs_info()

    pkgs = set()
    magics_re = None

    for cell in cells:
        if cell["cell_type"] != "code":
            continue

        # assuming we read the notebook with a json reader
        cell_source = "".join(cell["source"])
        if "import" not in cell_source:
            continue
        else:
            # quick hack to ignore jupyter magics
            if "%" in cell_source:
                if magics_re is None:
                    magics_re = re.compile(r"^( *)%{1,2}\w+ *", flags=re.MULTILINE)
                cell_source = magics_re.sub(r"\1", cell_source)

            for imp in cell_imports(cell_source):
                if imp in std_libs:  # type: ignore
                    continue
                if imp in pkgs_dists:  # type: ignore
                    pkgs.update(pkgs_dists[imp])  # type: ignore
                else:
                    pkgs.add(imp)

    if add_pkgs is not None:
        pkgs.update(add_pkgs)

    pkgs = {pkg: "" for pkg in pkgs}  # type: ignore
    if not pin_versions:
        return pkgs

    for pkg in pkgs:
        pkgs[pkg] = _get_version(pkg)  # type: ignore

    return pkgs


def _resolve(a, b):
    import packaging.version

    parse_version = packaging.version.parse

    if a == "":
        return b
    elif b == "":
        return a
    else:
        return a if parse_version(a) > parse_version(b) else b


def resolve_versions(notebooks_pkgs: list[dict]):
    """Harmonize packages' versions from lists of packages."""
    resolved = {}
    for pkgs in notebooks_pkgs:
        for pkg, ver in pkgs.items():
            if pkg not in resolved:
                resolved[pkg] = ver
            else:
                resolved[pkg] = _resolve(resolved[pkg], ver)

    return resolved
