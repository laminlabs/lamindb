from bionty import __version__ as bionty_v
from lndb_setup import __version__ as lndb_setup_v
from lnschema_core import __version__ as lnschema_core_v
from nbproject import __version__ as nbproject_v
from packaging import version

if version.parse(lndb_setup_v) != version.parse("0.27.0"):
    raise RuntimeError("lamindb needs lndb_setup==0.27.0")

if version.parse(lnschema_core_v) != version.parse("0.23.2"):
    raise RuntimeError("lamindb needs lnschema_core==0.23.2")

if version.parse(bionty_v) != version.parse("0.5.7"):
    raise RuntimeError("lamindb needs bionty==0.5.7")

if version.parse(nbproject_v) < version.parse("0.8.0"):
    raise RuntimeError("lamindb needs nbproject>=0.8.0")
