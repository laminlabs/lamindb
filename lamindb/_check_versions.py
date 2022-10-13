from bioreadout import __version__ as bioreadout_v
from lnbfx import __version__ as lnbfx_v
from lndb_hub import __version__ as lndb_hub_v
from lndb_setup import __version__ as lndb_setup_v
from lnschema_bionty import __version__ as lnschema_bionty_v
from lnschema_core import __version__ as lnschema_core_v
from lnschema_wetlab import __version__ as lnschema_wetlab_v
from nbproject import __version__ as nbproject_v
from packaging import version

if version.parse(lnschema_core_v) != version.parse("0.13.0"):
    raise RuntimeError("lamindb needs lnschema_core==0.13.0")

if version.parse(lnschema_wetlab_v) != version.parse("0.5.0"):
    raise RuntimeError("lamindb needs lnschema_wetlab==0.5.0")

if version.parse(lnschema_bionty_v) != version.parse("0.4.3"):
    raise RuntimeError("lamindb needs lnschema_bionty==0.4.3")

if version.parse(lndb_setup_v) != version.parse("0.12.1"):
    raise RuntimeError("lamindb needs lndb_setup==0.12.1")

if version.parse(lndb_hub_v) != version.parse("0.5.6"):
    raise RuntimeError("lamindb needs lndb_hub==0.5.6")

if version.parse(nbproject_v) < version.parse("0.7.0"):
    raise RuntimeError("lamindb needs nbproject>=0.7.0")

if version.parse(lnbfx_v) < version.parse("0.4.2"):
    raise RuntimeError("lamindb needs lnbfx>=0.4.2")

if version.parse(bioreadout_v) < version.parse("0.2.1"):
    raise RuntimeError("lamindb needs bioreadout>=0.2.1")
