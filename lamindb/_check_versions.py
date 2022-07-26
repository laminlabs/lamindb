from lndb_setup import __version__ as lndb_setup_version
from nbproject import __version__ as nbproject_version
from packaging import version

if version.parse(nbproject_version) < version.parse("0.4.3"):
    raise RuntimeError("lamindb needs nbproject >= 0.4.3")

if version.parse(lndb_setup_version) < version.parse("0.2.0"):
    raise RuntimeError("lamindb needs lndb_setup >= 0.2.0")
