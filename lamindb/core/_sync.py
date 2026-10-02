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
    """Sync objects from a source database to the default database.

    This function underlies the CLI command: `lamin io sync`.

    Guide: :doc:`transfer`

    Args:
        registry: Registry name, for example `artifact` or `record`.
        uid: UID of the object to sync.
        source_db: Source database slug, for example `laminlabs/lamindata`.
        depth: How many levels of records under a type to transfer.
            `0` transfers only the given object, plus related objects selected
            by `transfer`. A positive integer also transfers that many levels
            of records whose type chain starts at the object. Only `record`,
            `feature`, `schema`, `project`, `ulabel`, and `reference` accept
            `depth > 0`.
        transfer: `sqlrecord`, `notes`, or `annotations`.
            Omit it to use the registry default.

    One sync walks a single graph. Shared nodes are copied once. Link rows are
    replaced on the target, not copied by uid.

    .. code-block:: mermaid

       flowchart TD
         roots("Requested records and depth children") --> row("Copy the row once; fill Record and ULabel")
         row --> shared("Schema, features, dtype types: once per uid")
         row --> lookup("Annotation values: one uid lookup per registry")
         lookup --> present("Already on target: use that row")
         lookup --> stub("Missing Record or ULabel: stub")
         lookup --> once("Missing other registry: save once")
         row --> links("Replace link rows; do not copy them by uid")

    Most of the time, you will just the equivalent `.save()` on an object from another database::

        import lamindb as ln
        db = ln.DB("laminlabs/lamindata")
        # sync a record
        db.Record.get("gL3TbX2qZQmCwTAU").save()
        # sync a record with annotations
        db.Record.get("gL3TbX2qZQmCwTAU").save(transfer="annotations")
        # sync a record type and its data records at depth 1
        db.Record.get("gL3TbX2qZQmCwTAU").save(transfer="annotations", depth=1)
        # sync an artifact with annotations
        db.Artifact.get("gL3TbX2qZQmCwTAU").save(transfer="annotations")
        # sync a feature
        db.Feature.get("gL3TbX2qZQmCwTAU").save()
        # sync a schema
        db.Schema.get("gL3TbX2qZQmCwTAU").save()
        # sync a project
        db.Project.get("gL3TbX2qZQmCwTAU").save()
        # sync a ulabel
        db.ULabel.get("gL3TbX2qZQmCwTAU").save()

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
