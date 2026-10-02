from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ..models.sqlrecord import SQLRecord


def _registry_class_name(registry: str) -> str:
    if registry == "ulabel":
        return "ULabel"
    if not registry or not registry.replace("_", "").isalnum():
        raise ValueError(f"Unknown registry {registry!r}.")
    return "".join(part.capitalize() for part in registry.split("_"))


def sync(
    *,
    registry: str,
    uid: str,
    source_db: str,
    depth: int = 0,
    transfer: str | None = None,
) -> SQLRecord:
    """Sync one object from a source database into the current database.

    This function underlies `lamin io sync`.

    Guide: :doc:`transfer`

    Args:
        registry: Registry name, for example `artifact` or `record`.
        uid: UID of the object to sync.
        source_db: Source database slug, for example `laminlabs/lamindata`.
        depth: How many levels of the type tree to transfer. `0` transfers
            only this object, plus the related objects selected by `transfer`.
            Only `record`, `feature`, `schema`, `project`, `ulabel`, and
            `reference` accept `depth > 0`.
        transfer: `sqlrecord`, `notes`, or `annotations`.
            Omit it to use the registry default. Schema defaults to `annotations`.

    What is copied
    --------------

    `transfer` sets the boundary.

    .. code-block:: mermaid

       flowchart TD
         you("Object you sync") --> row("Its row and required foreign keys")
         row --> mode{"transfer"}
         mode --> bare("sqlrecord: stop after the row")
         mode --> notes("notes: also the latest readme")
         mode --> ann("annotations: also one step of links")
         ann --> feat("Features of this object")
         ann --> vals("Values linked from this object")
         vals --> stub("Record and ULabel: stub")
         vals --> other("Artifact and other registries: save with their annotations")
         you --> depth("depth, type tree only")
         depth --> kids("Direct records of this type, then depth - 1")

    `sqlrecord` copies the row. Foreign keys that the row needs are mapped by
    uid or created. `run` and `transform` are not copied from the source. They
    point at this transfer's run.

    `notes` also copies the latest readme.

    `annotations` also copies one step of links on this object: feature values,
    labels, and, for a schema, its members. It does not copy the annotations of
    those linked records.

    A linked `Record` or `ULabel` is a stub: uid, name, type, and creator. Its
    own features, labels, and readme stay on the source. Transfer that record
    itself, with `transfer="annotations"`, when you want them. A linked branch
    is the same kind of link. A stub is enough for a branch, because a branch
    has no annotations you are trying to keep. It is not enough for a record.

    A linked artifact, feature, schema, or other registry is saved with
    `transfer="annotations"`, so its own annotations come along. A data record
    whose type is not on the target yet is refused. Transfer that type first.

    `depth` only follows the type tree of `Record`, `Feature`, `Schema`,
    `Project`, `ULabel`, and `Reference`. `depth=1` adds the records whose type
    is the object you named. `depth=2` also adds the records typed by those.
    An artifact is never a depth child. An artifact is copied only when it is a
    foreign key or an annotation value of an object that is actually transferred.

    A record-frame is a record type. Its rows are data records of that type, so
    they are included only if you sync the type and pass `depth`. A row that
    merely appears as a feature value of something else is a stub: that sheet's
    other rows, and that row's own features, are not copied.

    .. code-block:: mermaid

       flowchart TD
         sheet("Sync the sheet type, depth=1") --> rows("Its rows are transferred")
         rows --> rowann("Each row keeps the transfer mode you passed")
         sample("Sync one sample") --> link("A feature points at a row of a sheet")
         link --> onerow("That one row is a stub")
         onerow --> notsheet("The rest of the sheet is not copied")

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
        db.Record.get("gL3TbX2qZQmCwTAU").save()
        db.Record.get("gL3TbX2qZQmCwTAU").save(transfer="annotations")
        db.Record.get("gL3TbX2qZQmCwTAU").save(transfer="annotations", depth=1)
        db.Artifact.get("gL3TbX2qZQmCwTAU").save(transfer="annotations")
    """
    from ..models.db import DB

    if type(depth) is not int or depth < 0:
        raise ValueError("depth must be an int >= 0.")
    model_name = _registry_class_name(registry)
    record = getattr(DB(source_db), model_name).get(uid)
    kwargs: dict[str, Any] = {"depth": depth}
    if transfer is not None:
        kwargs["transfer"] = transfer
    return record.save(**kwargs)
