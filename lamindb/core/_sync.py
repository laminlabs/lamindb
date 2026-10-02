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

    The `transfer` argument determines what is copied. If `"sqlrecord"`, only the object
    with its required foreign keys are copied. If `"notes"`, the object and its notes are copied.
    If `"annotations"`, the object and its annotations are copied, that is, the object's features,
    labels, and, for a schema, its members.

    .. code-block:: mermaid

       flowchart TD
         you("transfer = annotations") --> links("One step of links on this object")
         links --> feat("Features of this object")
         links --> vals("Values linked from this object")
         vals --> stub("Record and ULabel: stub")
         vals --> other("Artifact and other registries: save with their annotations")

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
        db.Record.get("gL3TbX2qZQmCwTAU").save(transfer="annotations")
        db.Artifact.get("gL3TbX2qZQmCwTAU").save(transfer="annotations")
    """
    if type(depth) is not int or depth < 0:
        raise ValueError("depth must be an int >= 0.")
    sqlrecord = registry.connect(source_db).get(uid)
    return sqlrecord.save(depth=depth, transfer=transfer)
