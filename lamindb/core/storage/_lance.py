from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Callable

    from lancedb.table import LanceTable
    from upath import UPath

    from lamindb.models.artifact import Artifact
    from lamindb.models.run import Run


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


class LanceAccessor:
    """Read and revise a `.lance` artifact.

    Reached through :attr:`~lamindb.Artifact.lance`. :meth:`open` returns the
    LanceDB table at the version pinned on that artifact. :meth:`add` and
    :meth:`add_columns` write to the table and save a new artifact version.
    The artifact you called them on stays pinned to the Lance version from
    before the write.
    """

    def __init__(self, artifact: Artifact) -> None:
        if artifact.suffix != ".lance":
            raise ValueError("artifact.lance expects a .lance artifact.")
        self._artifact = artifact

    def __repr__(self) -> str:
        artifact = self._artifact
        return f"LanceAccessor({artifact.key or artifact.uid})"

    def open(self, *, version: int | None = None) -> LanceTable:
        """Open this artifact as one LanceDB table.

        Returns a :class:`lancedb.table.LanceTable`. When this artifact version
        has a pinned Lance version and that version is not the directory's
        latest, the table is checked out at the pin. Pass `version` to check
        out a specific Lance version instead.

        Writes on the returned table change the Lance directory and leave this
        artifact's pin unchanged. :meth:`add` and :meth:`add_columns` write and
        save a new artifact version.

        On S3, AWS credentials are read from the artifact path when the table
        is opened. A refreshable session is snapshotted at that moment.

        Args:
            version: Lance version to check out. Defaults to the version pinned
                when this artifact version was written.
        """
        artifact = self._artifact
        path = getattr(artifact, "_local_filepath", None)
        if path is None or not path.exists():
            path = artifact.path
        table = open_lance_table(path)
        pinned = lance_version(artifact)
        if version is not None:
            table.checkout(version)
        elif pinned is not None and pinned != int(table.version):
            table.checkout(pinned)
        return table

    def add(self, data: Any, *, run: Run | None = None) -> Artifact:
        """Append rows and save a new artifact version.

        Calls `LanceTable.add(data, mode="append")` on this artifact's
        directory, then registers a new artifact version of that same directory.
        The previous version stays pinned to the Lance version from before the
        append. `lance.open()` checks out that version.

        The directory is shared across these artifact versions, as with a
        `.tiledbsoma` store. `from_dataframe()` still writes a new directory.
        Writes stay on Lance's default branch. Coupling LaminDB branches to
        Lance branches is not supported and not needed for isolation: a version
        saved on another LaminDB branch is not visible on main.

        Args:
            data: Rows to append. Any object accepted by `LanceTable.add`.
            run: The run that creates the new artifact version.

        Examples:

            Append rows and keep the previous Lance version readable::

                revised = artifact.lance.add(new_rows)
                artifact.lance.open().count_rows()
                revised.lance.open().count_rows()
        """

        def append_rows(table: LanceTable) -> None:
            table.add(data, mode="append")

        return self._revise(append_rows, run=run)

    def add_columns(
        self,
        transforms: Any = None,
        *,
        computed: dict[str, str] | None = None,
        run: Run | None = None,
    ) -> Artifact:
        """Add columns and save a new artifact version.

        Calls `LanceTable.add_columns` on this artifact's directory, then
        registers a new artifact version of that same directory. The previous
        version stays pinned to the Lance version from before the write, so
        opening it still returns the previous schema.

        Args:
            transforms: Column definitions accepted by `LanceTable.add_columns`.
                A dict maps a new column name to a SQL expression evaluated on
                each existing row. A pyarrow field, list of fields, or schema
                adds columns filled with nulls.
            computed: SQL expressions for maintained computed columns.
            run: The run that creates the new artifact version.

        Examples:

            Add a column computed from existing rows::

                revised = artifact.lance.add_columns({"score": "n_genes * 2"})
        """

        def append_columns(table: LanceTable) -> None:
            table.add_columns(transforms, computed=computed)

        return self._revise(append_columns, run=run)

    def _revise(
        self, write: Callable[[LanceTable], None], *, run: Run | None
    ) -> Artifact:
        """Apply `write` to the table and save the result as a new version."""
        artifact = self._artifact
        if artifact._state.adding:
            raise ValueError("Save the artifact before revising it.")
        if not artifact.is_latest:
            raise ValueError("Write to the latest version of the artifact.")

        from lamindb.models._lineage import get_run
        from lamindb.models.artifact import Artifact

        # Capture the path before flipping overwrite_versions. A version with
        # the flag unset is stored under its full uid; the shared directory
        # uses the 16-character stem. `.lance` defaults to the stem.
        store_path = artifact.path
        table = open_lance_table(store_path)
        if lance_version(artifact) is None:
            set_lance_version(artifact, int(table.version))
        write(table)
        n_observations = table.count_rows()
        new_lance_version = int(table.version)
        del table

        stem_path = store_path.with_name(f"{artifact.stem_uid}{artifact.suffix}")
        if store_path != stem_path:
            store_path.rename(stem_path)
        if not artifact._overwrite_versions:
            artifact._overwrite_versions = True
        artifact.save()

        revised = Artifact(  # type: ignore
            stem_path,
            key=artifact.key,
            description=artifact.description,
            kind=artifact.kind,
            run=get_run(run),
            revises=artifact,
            otype=artifact.otype,
            overwrite_versions=True,
            _is_internal_call=True,
        )
        revised.n_observations = n_observations
        set_lance_version(revised, new_lance_version)
        return revised.save()
