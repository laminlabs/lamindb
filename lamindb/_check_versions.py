from bioreadout import __version__ as bioreadout_v
from lnbfx import __version__ as lnbfx_v
from lndb_hub import __version__ as lndb_hub_v
from lndb_setup import __version__ as lndb_setup_v
from lnschema_core import __version__ as lnschema_core_v
from lnschema_wetlab import __version__ as lnschema_wetlab_v
from nbproject import __version__ as nbproject_v
from packaging import version

if version.parse(lnschema_core_v) != version.parse("0.5.1"):
    raise RuntimeError("lamindb needs lnschema_core==0.5.1")

if version.parse(lnschema_wetlab_v) != version.parse("0.3.7"):
    raise RuntimeError("lamindb needs lnschema_wetlab==0.3.7")

if version.parse(lndb_setup_v) != version.parse("0.6.3"):
    raise RuntimeError("lamindb needs lndb_setup==0.6.3")

if version.parse(lndb_hub_v) != version.parse("0.5.6"):
    raise RuntimeError("lamindb needs lndb_hub==0.5.6")

if version.parse(nbproject_v) < version.parse("0.5.0"):
    raise RuntimeError("lamindb needs nbproject>=0.5.0")

if version.parse(lnbfx_v) < version.parse("0.3.4"):
    raise RuntimeError("lamindb needs lnbfx>=0.3.4")

if version.parse(bioreadout_v) != version.parse("0.1.0"):
    raise RuntimeError("lamindb needs bioreadout==0.1.0")
