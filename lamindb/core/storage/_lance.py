from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from lancedb.table import LanceTable
    from upath import UPath


# Lance version pinned on the artifact that wrote it. Stored in `Artifact._aux` as "lv".
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


def open_lance_table(path: UPath) -> LanceTable:
    """Open one `.lance` directory as a LanceDB table.

    Connects to the parent of this directory and opens the table by name.
    Sibling directories are not part of the table.
    """
    lancedb = import_lancedb()
    database = lancedb.connect(path.parent.as_posix())
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
