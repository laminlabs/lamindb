from typing import Optional

from lamin_logger import logger
from lnschema_core import DObject
from sqlalchemy.orm.session import object_session

from lamindb._context import context

from ._settings import settings
from .dev._core import filepath_from_dobject
from .dev.db import add
from .dev.file import load_to_memory


# this is exposed to the user as DObject.load
def load(dobject: DObject, stream: bool = False, is_run_input: Optional[bool] = None):
    if stream and dobject.suffix not in (".h5ad", ".zarr"):
        logger.warning(f"Ignoring stream option for a {dobject.suffix} object.")
    if is_run_input is None:
        track_run_input = settings.track_run_inputs_upon_load
    else:
        track_run_input = is_run_input
    if track_run_input:
        if object_session(dobject) is None:
            raise ValueError("Need to load with session open to track as input.")
        if context.run is None:
            raise ValueError(
                "No global run context set. Call ln.context.track() or pass input run"
                " directly."
            )
        else:
            dobject.targets.append(context.run)
            add(dobject)  # need to commit here
    # track_usage(dobject.id, "load")
    return load_to_memory(filepath_from_dobject(dobject), stream=stream)
