from typing import Optional

from lamin_logger import logger
from lndb_storage import load_to_memory
from lnschema_core import File
from sqlalchemy.orm.session import object_session

from lamindb._context import context

from ._settings import settings
from .dev._core import filepath_from_file


# this is exposed to the user as File.load
def load(file: File, stream: bool = False, is_run_input: Optional[bool] = None):
    if stream and file.suffix not in (".h5ad", ".zarr"):
        logger.warning(f"Ignoring stream option for a {file.suffix} object.")
    if is_run_input is None:
        track_run_input = settings.track_run_inputs_upon_load
    else:
        track_run_input = is_run_input
    if track_run_input:
        if object_session(file) is None:
            raise ValueError("Need to load with session open to track as input.")
        if context.run is None:
            raise ValueError(
                "No global run context set. Call ln.context.track() or pass input run"
                " directly."
            )
        else:
            file.targets.append(context.run)
            session = object_session(file)
            session.add(file)
            session.commit()
    return load_to_memory(filepath_from_file(file), stream=stream)
