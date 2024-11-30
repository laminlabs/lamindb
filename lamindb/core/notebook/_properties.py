from __future__ import annotations

from typing import TYPE_CHECKING

from lamin_utils import logger

if TYPE_CHECKING:
    from nbformat import NotebookNode


def get_title(nb: NotebookNode) -> str | None:
    """Get title of the notebook."""
    # loop through all cells
    for cell in nb.cells:
        # only consider markdown
        if cell["cell_type"] == "markdown":
            # grab source
            text = "".join(cell["source"])
            # loop through lines
            for line in text.split("\n"):
                # if finding a level-1 heading, consider it a title
                if line.startswith("# "):
                    title = line.lstrip("#").strip(" .").strip("\n")
                    return title
    return None


def check_consecutiveness(nb: NotebookNode, calling_statement: str = None) -> bool:
    """Check whether code cells have been executed consecutively.

    Needs to be called in the last code cell of a notebook.
    Otherwise raises `RuntimeError`.

    Returns cell transitions that violate execution at increments of 1 as a list
    of tuples.

    Args:
        nb: Notebook content.
        calling_statement: The statement that calls this function.
    """
    cells = nb.cells

    violations = []
    prev = 0

    ccount = 0  # need to initialize because notebook might note have code cells
    # and below, we check if ccount is None
    for cell in cells:
        cell_source = "".join(cell["source"])
        if cell["cell_type"] != "code" or cell_source == "":
            continue

        if calling_statement is not None and calling_statement in cell_source:
            continue

        ccount = cell["execution_count"]
        if ccount is None or prev is None or ccount - prev != 1:
            violations.append((prev, ccount))

        prev = ccount

    # ignore the very last code cell of the notebook
    # `check_consecutiveness` is being run during publish if `last_cell`` is True
    # hence, that cell has ccount is None
    if ccount is None:
        violations.pop()

    any_violations = len(violations) > 0
    if any_violations:
        logger.warning(f"cells {violations} were not run consecutively")
    else:
        logger.success("cell execution numbers increase consecutively")

    return not any_violations
