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

    This function underlies `lamin io sync`.

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

    **Example:** For an artifact, the following relationships are foreign keys, which are copied even
    when `transfer="sqlrecord"`.

    .. code-block:: mermaid

       flowchart TD
         artifact("Artifact") -->|storage| storage("Storage")
         artifact -->|branch| branch("Branch")
         artifact -->|created_by| user("User")
         artifact -->|created_on| created_on("Branch")
         artifact -->|run| run("Run")
         artifact -->|schema| schema("Schema")
         artifact -->|space| space("Space")

    A linked `Record` or `ULabel` is a stub: uid, name, type, and creator. Its
    own features, labels, and readme stay on the source. Transfer that record
    itself when you want them. A linked branch is a stub, and a stub is enough
    for a branch.

    A linked artifact, feature, schema, or other registry is saved with
    `transfer="annotations"`, so its own annotations come along. A data record
    whose type is not on the target yet is refused. Transfer that type first.

    `sqlrecord` stops after the row. `notes` also copies the latest readme.
    `run` and `transform` point at this transfer's run, not the source run.

    Running it again
    ----------------

    A transfer is safe to repeat. Uids already on the target are reused. A
    `Record` or `ULabel` that arrived earlier as a stub is filled in when you
    transfer that object itself. Link rows are replaced, not duplicated. So a
    first run with `transfer="sqlrecord"` and a second run with
    `transfer="annotations"` completes the annotations.

    Every row this transfer writes has `.run` set to a run of the transform
    `__lamindb_transfer__/{source instance uid}`. To undo it, find that run and
    delete the objects whose `.run` is that run. Objects that were already on
    the target and only got mapped are not part of that run.

    Most of the time, call `.save()` on the object from the other database::

        import lamindb as ln

        db = ln.DB("laminlabs/lamindata")
        db.Artifact.get(key="example_datasets/mini_immuno/dataset1.h5ad").save(
            transfer="annotations"
        )
        db.Record.get("gL3TbX2qZQmCwTAU").save(transfer="annotations")
    """
    if type(depth) is not int or depth < 0:
        raise ValueError("depth must be an int >= 0.")
    sqlrecord = registry.connect(source_db).get(uid)
    return sqlrecord.save(depth=depth, transfer=transfer)
