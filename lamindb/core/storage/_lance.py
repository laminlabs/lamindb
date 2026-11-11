from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from lancedb.table import LanceTable
    from upath import UPath


def open_lance_table(path: UPath, **kwargs: Any) -> LanceTable:
    """Open one `.lance` directory as a LanceDB table.

    Connects to the parent of this directory and opens the table by name.
    Sibling directories are not part of the table.
    """
    try:
        import lancedb
    except ImportError as exc:
        raise ImportError(
            "Please install lancedb: pip install 'lamindb[lance]'"
        ) from exc

    connect_kwargs: dict[str, Any] = {}
    if "storage_options" in kwargs:
        connect_kwargs["storage_options"] = kwargs.pop("storage_options")
    else:
        storage_options = lance_storage_options(path)
        if storage_options:
            connect_kwargs["storage_options"] = storage_options
    database = lancedb.connect(path.parent.as_posix(), **connect_kwargs)
    return database.open_table(path.stem, **kwargs)


def lance_storage_options(path: UPath) -> dict[str, str] | None:
    """Map a cloud filesystem's credentials onto LanceDB storage options."""
    protocol = getattr(path, "protocol", None)
    if protocol in {None, "file", "local", ""}:
        return None

    from ._polars_lazy_df import _polars_options

    polars_options = _polars_options(path)
    storage = dict(polars_options.get("storage_options") or {})
    provider = polars_options.get("credential_provider")
    if provider is not None:
        credentials, _expiry = provider()
        storage.update(credentials)

    options: dict[str, str] = {}
    for source, target in (
        ("aws_access_key_id", "aws_access_key_id"),
        ("aws_secret_access_key", "aws_secret_access_key"),
        ("aws_session_token", "aws_session_token"),
        ("aws_region", "aws_region"),
    ):
        value = storage.get(source)
        if value:
            options[target] = value
    endpoint = storage.get("aws_endpoint_url")
    if endpoint:
        options["aws_endpoint"] = endpoint
        if storage.get("aws_allow_http"):
            options["allow_http"] = "true"
    if storage.get("aws_skip_signature"):
        options["aws_skip_signature"] = "true"
    return options or None
