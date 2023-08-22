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
    if run is not None:
        kwargs["transform"] = run.transform
