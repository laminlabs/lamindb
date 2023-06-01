import os

from lamindb_setup import __version__ as lndb_v
from lndb_storage import __version__ as lndb_storage_v
from lnschema_core import __version__ as lnschema_core_v
from nbproject import __version__ as nbproject_v
from packaging import version

if os.getenv("GITHUB_ACTIONS") is None:
    # Lamin PINNED packages

    if version.parse(lnschema_core_v) != version.parse("0.34.0"):
        raise RuntimeError("lamindb needs lnschema_core==0.34.0")

    if version.parse(lndb_storage_v) != version.parse("0.3.2"):
        raise RuntimeError("lamindb needs lndb_storage==0.3.2")

    if version.parse(lndb_v) < version.parse("0.45.0"):
        raise RuntimeError("Upgrade lndb! pip install lndb==0.45.0")

    # Lamin GREATEREQ packages

    if version.parse(nbproject_v) < version.parse("0.8.5"):
        raise RuntimeError("lamindb needs nbproject>=0.8.5")
