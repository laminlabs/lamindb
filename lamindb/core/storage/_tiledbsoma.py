from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from lamindb_setup.core._settings_storage import get_storage_region
from lnschema_core import Artifact

if TYPE_CHECKING:
    from tiledbsoma import Collection as SOMACollection
    from tiledbsoma import Experiment as SOMAExperiment
    from upath import UPath


def _tiledb_config_s3(filepath: UPath) -> dict:
    region = get_storage_region(filepath)
    tiledb_config = {"vfs.s3.region": region}
    storage_options = filepath.storage_options
    if "key" in storage_options:
        tiledb_config["vfs.s3.aws_access_key_id"] = storage_options["key"]
    if "secret" in storage_options:
        tiledb_config["vfs.s3.aws_secret_access_key"] = storage_options["secret"]
    if "token" in storage_options:
        tiledb_config["vfs.s3.aws_session_token"] = storage_options["token"]

    return tiledb_config


def _open_tiledbsoma(
    filepath: UPath, mode: Literal["r", "w"] = "r"
) -> SOMACollection | SOMAExperiment:
    try:
        import tiledbsoma as soma
    except ImportError as e:
        raise ImportError("Please install tiledbsoma: pip install tiledbsoma") from e

    filepath_str = filepath.as_posix()
    if filepath.protocol == "s3":
        ctx = soma.SOMATileDBContext(tiledb_config=_tiledb_config_s3(filepath))
        # this is a strange bug
        # for some reason iterdir futher gives incorrect results
        # if cache is not invalidated
        # instead of obs and ms it gives ms and ms in the list of names
        filepath.fs.invalidate_cache()
    else:
        ctx = None

    soma_objects = [obj.name for obj in filepath.iterdir()]
    if "obs" in soma_objects and "ms" in soma_objects:
        SOMAType = soma.Experiment
    else:
        SOMAType = soma.Collection
    return SOMAType.open(filepath_str, mode=mode, context=ctx)
