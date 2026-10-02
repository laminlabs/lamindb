from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..base.types import TransferMode
    from ..models.sqlrecord import Registry, SQLRecord


def sync(
    *,
    registry: Registry,
    uid: str,
    source_db: str,
    transfer: TransferMode | None = None,
    depth: int = 0,
) -> SQLRecord:
    """Sync one object from a source database into the current database.

    This function underlies `lamin io sync` and wraps the lower-level `.save()` API::

        import lamindb as ln

        db = ln.DB("laminlabs/lamindata")
        db.Record.get("gL3TbX2qZQmCwTAU").save(transfer="annotations")

    Guide: :doc:`transfer`

    Args:
        registry: Registry class, for example `ln.Artifact` or `ln.Record`.
        uid: UID of the object to sync.
        source_db: Source database slug, for example `laminlabs/lamindata`.
        transfer: A :class:`~lamindb.base.types.TransferMode`.
            Omit it to use the registry default. Schema defaults to `annotations`.
        depth: How many levels of the type tree to transfer. `0` transfers
            only this object, plus the related objects selected by `transfer`.
            Only `record`, `feature`, `schema`, `project`, `ulabel`, and
            `reference` accept `depth > 0`.

    Returns:
        The saved `SQLRecord` object on the current database.

    Traversing relationships
    ------------------------

    The `transfer` argument determines which related objects are transferred.
    If `"sqlrecord"`, only the object with its required foreign keys are copied.
    If `"notes"`, the object and its notes are copied.
    If `"annotations"`, the object and its annotations are copied, that is, the object's features,
    labels, and, for a schema, its members.

    **Example:** If you pass `transfer="sqlrecord"` upon transferring an artifact, the following foreign keys are transferred:

    .. code-block:: mermaid

       flowchart TD
         artifact("Artifact") -->|storage| storage("Storage")
         artifact -->|branch| branch("Branch")
         artifact -->|created_by| user("User")
         artifact -->|created_on| created_on("Branch")
         artifact -->|run| run("Run")
         artifact -->|schema| schema("Schema")
         artifact -->|space| space("Space")

    If you pass `transfer="annotations"` upon transferring an artifact, its many-to-many relationships are transferred.
    Those contain all label & feature annotations but also inferred schemas via `.schemas`.
    Some of the many-to-many relationships of an artifact are shown below:

    .. code-block:: mermaid

       flowchart TD
         artifact("Artifact") -->|ulabels| ulabels("ULabel")
         artifact -->|records| records("Record")
         artifact -->|projects| projects("Project")
         artifact -->|users| users("User")
         artifact -->|artifacts| linked("Artifact")
         artifact -->|schemas| schemas("Schema")
         artifact -->|json_values| json_values("JsonValue")

    The related objects themselves are **transferred** without their own annotations to avoid an infinite recursion.
    You have to transfer the related object itself if you want to transfer it with its own annotations.

    Re-syncing
    ----------

    A sync operation is safe to repeat. UIDs already on the target database are mapped. A
    `Record` that arrived earlier as a stub is filled in when you
    transfer that object itself. Links are replaced, not duplicated. So a
    first run with `transfer="sqlrecord"` and a second run with
    `transfer="annotations"` completes annotations.

    Every row this transfer writes has `.run` set to a run of the transform
    `__lamindb_transfer__/{source_database_uid}`. To undo it, find that run and
    delete the objects whose `.run` is that run. Objects that were already on
    the target and only got mapped are not part of that run.
    """
    if type(depth) is not int or depth < 0:
        raise ValueError("depth must be an int >= 0.")
    sqlrecord = registry.connect(source_db).get(uid)
    return sqlrecord.save(depth=depth, transfer=transfer)
