import os

from lamindb_setup import __version__ as lndb_v
from lnschema_core import __version__ as lnschema_core_v
from packaging import version

if os.getenv("GITHUB_ACTIONS") is None:
    # Lamin PINNED packages

    if version.parse(lnschema_core_v) != version.parse("0.35a4"):
        raise RuntimeError("lamindb needs lnschema_core==0.35a4")

    if version.parse(lndb_v) < version.parse("0.46a3"):
        raise RuntimeError("Upgrade lndb! pip install lndb==0.46a3")

    # Lamin GREATEREQ packages
    try:
        from nbproject import __version__ as nbproject_v

        if version.parse(nbproject_v) < version.parse("0.8.7"):
            raise RuntimeError("lamindb needs nbproject>=0.8.7")
    except ImportError:
        pass
