from typing import Any, Dict, Optional

from lamin_utils import logger
from lnschema_core import Run

from lamindb.dev import run_context


def get_run(run: Optional[Run]) -> Optional[Run]:
    if run is None:
        run = run_context.run
        if run is None:
            logger.warning(
                "no run & transform get linked, consider passing a `run` or calling"
                " ln.track()"
            )
    return run


def add_transform_to_kwargs(kwargs: Dict[str, Any], run: Run):
    # transform cannot be directly passed, just via run
    # it's directly stored in the file table to avoid another join
    # mediate by the run table
    if run is not None:
        if run.transform_id is not None:
            kwargs["transform_id"] = run.transform_id
        else:
            # accessing the relationship should always be possible if
            # the above if clause was false as then, we should have a fresh
            # Transform object that is not queried from the DB
            assert run.transform is not None
            kwargs["transform"] = run.transform
