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
    from lamin_cli._save import save_run_context_core

    if is_run_from_ipython:
        # notebooks
        from nbproject.dev import read_notebook
        from nbproject.dev._check_last_cell import check_last_cell

        if not i_saved_the_notebook and not ln_setup._TESTING:
            logger.error(
                "Save the notebook, pass `i_saved_the_notebook=True`, and re-run this cell."
            )
            return None
        notebook_content = read_notebook(run_context.path)  # type: ignore
        if not check_last_cell(notebook_content, "i_saved_the_notebook"):
            raise CallFinishInLastCell(
                "Can only finish() from the last code cell of the notebook."
            )
        save_run_context_core(
            run=run_context.run,
            transform=run_context.transform,
            filepath=run_context.path,
            notebook_content=notebook_content,
            finished_at=True,
        )
    else:
        # scripts
        run_context.run.finished_at = datetime.now(timezone.utc)  # update run time
        run_context.run.save()
