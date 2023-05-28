from lamin_logger import logger
from lndb import __version__ as lndb_v
from lndb_storage import __version__ as lndb_storage_v
from lnschema_core import __version__ as lnschema_core_v
from nbproject import __version__ as nbproject_v
from packaging import version

# Lamin PINNED packages

if version.parse(lnschema_core_v) != version.parse("0.34a2"):
    raise RuntimeError("lamindb needs lnschema_core==0.34a2")

if version.parse(lndb_storage_v) != version.parse("0.3.2"):
    raise RuntimeError("lamindb needs lndb_storage==0.3.2")

if version.parse(lndb_v) < version.parse("0.45a4"):
    raise RuntimeError("Upgrade lndb! pip install lndb==0.45a4")

# Lamin GREATEREQ packages

if version.parse(nbproject_v) < version.parse("0.8.5"):
    raise RuntimeError("lamindb needs nbproject>=0.8.5")

# ensure that the lamin package is not installed
try:
    import lamin  # noqa

    logger.warning(
        "Please,\n"
        " - replace `import lamin` with `import lamindb.setup as lnsetup`\n"
        " - run `pip uninstall lamin`\n"
        "lamindb.setup now has all of the lamin functionality\n"
        "The lamindb API and lamin API will be integrated soon!\n"
        "The CLI remains as is!"
    )
except ImportError:
    pass
