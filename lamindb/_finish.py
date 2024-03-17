from datetime import datetime, timezone

import lamindb_setup as ln_setup
from lamin_utils import logger

from .core._run_context import is_run_from_ipython, run_context


class CallFinishInLastCell(SystemExit):
    pass


def finish(i_saved_the_notebook: bool = False):
    """Mark the tracked run as finished.

    Save the run report to your default storage location.
    """
    from lamin_cli._save import save

    if is_run_from_ipython:
        from nbproject.dev import read_notebook
        from nbproject.dev._check_last_cell import check_last_cell

        if not i_saved_the_notebook and not ln_setup._TESTING:
            logger.error(
                "Save the notebook, pass `i_saved_the_notebook=True`, and re-run this cell."
            )
            return None
        nb = read_notebook(run_context.path)  # type: ignore
        if not check_last_cell(nb, "i_saved_the_notebook"):
            raise CallFinishInLastCell(
                "Can only finish() from the last code cell of the notebook."
            )
        # scripts are already saved during `ln.track()`
        # TODO: make this more symmetric
        save(run_context.path)

    run_context.run.finished_at = datetime.now(timezone.utc)  # update run time
    run_context.run.save()
