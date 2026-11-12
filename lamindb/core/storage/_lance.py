from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from lancedb.table import LanceTable
    from upath import UPath


# Lance version pinned on the artifact that wrote it. Stored in `Artifact._aux` as "lv".
# This is a version on Lance's default branch. Coupling LaminDB branches to Lance
# branches is not supported and not needed for isolation: an artifact version on
# LaminDB main points at the latest Lance version visible on main, and a version
# saved on another LaminDB branch is not visible there.
LANCE_VERSION_AUX_KEY = "lv"


def lance_version(artifact: Any) -> int | None:
    """Return the Lance version pinned on an artifact, if one was recorded."""
    aux = getattr(artifact, "_aux", None) or {}
    value = aux.get(LANCE_VERSION_AUX_KEY)
    if value is None:
        return None
    return int(value)


def set_lance_version(artifact: Any, version: int) -> None:
    """Pin a Lance version on an artifact."""
    aux = dict(getattr(artifact, "_aux", None) or {})
    aux[LANCE_VERSION_AUX_KEY] = int(version)
    artifact._aux = aux


def record_current_lance_version(artifact: Any) -> None:
    """Pin whatever Lance version the artifact's directory is currently at."""
    if lance_version(artifact) is not None:
        return
    path = getattr(artifact, "_local_filepath", None) or artifact.path
    set_lance_version(artifact, int(open_lance_table(path).version))


def import_lancedb() -> Any:
    """Import `lancedb` or raise the install hint for the `lance` extra."""
    try:
        import lancedb
    except ImportError as exc:
        raise ImportError(
            "Please install lancedb: pip install 'lamindb[lance]'"
        ) from exc
    return lancedb


def lance_storage_options(path: UPath) -> dict[str, str] | None:
    """Map an S3 filesystem's credentials onto LanceDB storage options.

    Uses the same filesystem credentials as the Polars reader. Lance reads
    these options once, when the table is opened, so a refreshable session is
    snapshotted at that moment.
    """
    if getattr(path, "protocol", None) != "s3":
        return None

    from ._polars_lazy_df import _polars_options

    polars_options = _polars_options(path)
    storage = dict(polars_options.get("storage_options") or {})
    provider = polars_options.get("credential_provider")
    if provider is not None:
        credentials, _expiry = provider()
        storage.update(credentials)

    options: dict[str, str] = {}
    for source in (
        "aws_access_key_id",
        "aws_secret_access_key",
        "aws_session_token",
        "aws_region",
    ):
        value = storage.get(source)
        if value:
            options[source] = str(value)
    endpoint = storage.get("aws_endpoint_url")
    if endpoint:
        options["aws_endpoint"] = str(endpoint)
        # A custom endpoint is path-style, matching the Polars reader.
        options["aws_virtual_hosted_style_request"] = "false"
        if storage.get("aws_allow_http"):
            options["allow_http"] = "true"
    if storage.get("aws_skip_signature"):
        options["skip_signature"] = "true"
    return options or None


def open_lance_table(path: UPath) -> LanceTable:
    """Open one `.lance` directory as a LanceDB table.

    Connects to the parent of this directory and opens the table by name.
    Sibling directories are not part of the table. An S3 path receives the
    filesystem's AWS credentials.
    """
    lancedb = import_lancedb()
    connect_kwargs: dict[str, Any] = {}
    storage_options = lance_storage_options(path)
    if storage_options:
        connect_kwargs["storage_options"] = storage_options
    database = lancedb.connect(path.parent.as_posix(), **connect_kwargs)
    return database.open_table(path.stem)


def write_lance_table(data: Any, filepath: Path | str, **kwargs: Any) -> None:
    """Write one LanceDB table as a `.lance` directory."""
    lancedb = import_lancedb()
    path = Path(filepath)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_kwargs = dict(kwargs)
    # A new artifact is a new table. Appends belong to a later version.
    write_kwargs["mode"] = "create"
    database = lancedb.connect(path.parent.as_posix())
    database.create_table(path.stem, data, **write_kwargs)
