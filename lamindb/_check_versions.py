from lndb_cli import __version__ as lndb_cli_version
from nbproject import __version__ as nbproject_version
from packaging import version

if version.parse(nbproject_version) < version.parse("0.4.3"):
    raise RuntimeError("lamindb needs nbproject >= 0.4.3")

if version.parse(lndb_cli_version) < version.parse("0.1.1"):
    raise RuntimeError("lamindb needs lndb_cli >= 0.1.1")
