from datetime import datetime, timezone

from lamin_utils import logger
from nbproject.dev import read_notebook
from nbproject.dev._check_last_cell import check_last_cell

from .core._run_context import is_run_from_ipython, run_context


def finish(i_saved_the_notebook: bool = False):
    """Mark the transform run as finished, upload notebook."""
    from lamin_cli._save import save

    if is_run_from_ipython:
        if not i_saved_the_notebook:
            logger.error(
                "Save the notebook in your editor before finishing and pass `i_saved_the_notebook=True`"
            )
            return None
        nb = read_notebook(run_context.path)  # type: ignore
        if not check_last_cell(nb, "i_saved_the_notebook"):
            logger.error("Can only finish() from the last code cell of the notebook.")
            return None

    save(run_context.path)
    run_context.run.finished_at = datetime.now(timezone.utc)  # update run time
    run_context.run.save()
