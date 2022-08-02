from lndb_schema_core import __version__ as lndb_schema_core_v
from lndb_setup import __version__ as lndb_setup_v
from nbproject import __version__ as nbproject_v
from packaging import version

if version.parse(nbproject_v) < version.parse("0.4.3"):
    raise RuntimeError("lamindb needs nbproject >= 0.4.3")

if version.parse(lndb_setup_v) < version.parse("0.4.2"):
    raise RuntimeError("lamindb needs lndb_setup >= 0.4.2")

if version.parse(lndb_schema_core_v) < version.parse("0.2.0"):
    raise RuntimeError("lamindb needs lndb_setup >= 0.2.0")
