from lndb_hub import __version__ as lndb_hub_v
from lndb_setup import __version__ as lndb_setup_v
from lnschema_core import __version__ as lnschema_core_v
from nbproject import __version__ as nbproject_v
from packaging import version

if version.parse(nbproject_v) < version.parse("0.4.3"):
    raise RuntimeError("lamindb needs nbproject >= 0.4.3")

if version.parse(lndb_setup_v) < version.parse("0.5.4"):
    raise RuntimeError("lamindb needs lndb_setup >= 0.5.4")

if version.parse(lndb_hub_v) < version.parse("0.4.0"):
    raise RuntimeError("lamindb needs lndb_hub >= 0.4.0")

if version.parse(lnschema_core_v) < version.parse("0.4.0"):
    raise RuntimeError("lamindb needs lnschema_core >= 0.4.0")
