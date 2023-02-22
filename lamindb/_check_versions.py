from bionty import __version__ as bionty_v
from lndb import __version__ as lndb_v
from lnschema_core import __version__ as lnschema_core_v
from nbproject import __version__ as nbproject_v
from packaging import version

if version.parse(lndb_v) != version.parse("0.35.2"):
    raise RuntimeError("Upgrade lndb! pip install lndb==0.35.2")

if version.parse(lnschema_core_v) != version.parse("0.28.0"):
    raise RuntimeError("lamindb needs lnschema_core==0.28.0")

if version.parse(bionty_v) != version.parse("0.7.0"):
    raise RuntimeError("lamindb needs bionty==0.7.0")

if version.parse(nbproject_v) < version.parse("0.8.2"):
    raise RuntimeError("lamindb needs nbproject>=0.8.2")
