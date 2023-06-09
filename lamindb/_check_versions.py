import os

from lamin_logger import logger
from lamindb_setup import __version__ as lamindb_setup_v
from lnschema_core import __version__ as lnschema_core_v
from packaging import version

if os.getenv("GITHUB_ACTIONS") is None:
    # Lamin PINNED packages

    if version.parse(lnschema_core_v) != version.parse("0.35.2"):
        logger.warning("lamindb needs lnschema_core==0.35.2")

    if version.parse(lamindb_setup_v) < version.parse("0.47.1"):
        logger.warning("Upgrade lamindb_setup! pip install lamindb_setup==0.47.1")

    # Lamin GREATEREQ packages
    try:
        from nbproject import __version__ as nbproject_v

        if version.parse(nbproject_v) < version.parse("0.8.7"):
            logger.warning("lamindb needs nbproject>=0.8.7")
    except ImportError:
        pass
