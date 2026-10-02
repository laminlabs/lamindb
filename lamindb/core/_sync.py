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
    registry: str,
    uids: str | list[str],
    *,
    source: str,
    depth: int = 0,
    transfer: str | None = None,
) -> list[SQLRecord]:
    """Sync SQLRecord objects from a source database into the default database.

    This is a high-level function used in the CLI: `lamin io sync`.

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

    Most of the time, you will just `.save()` on an object from another database::

        import lamindb as ln
        db = ln.DB("laminlabs/lamindata")
        record = db.Record.get(uid="gL3TbX2qZQmCwTAU")
        record.save(transfer="sqlrecord")

    Guide: :doc:`transfer`

    Args:
        registry: Registry name, for example `artifact` or `record`.
        uids: One uid or several uids on the source database.
        source: Source instance slug, for example `laminlabs/lamindata`.
        depth: How many levels of records under a type to transfer.
            `0` transfers only the given objects, plus related objects selected
            by `transfer`. A positive integer also transfers that many levels
            of records whose type chain starts at each object. Only `record`,
            `feature`, `schema`, `project`, `ulabel`, and `reference` accept
            `depth > 0`.
        transfer: `sqlrecord`, `notes`, or `annotations`.
            Omit it to use the registry default.
    """
    from ..models.db import DB

    if type(depth) is not int or depth < 0:
        raise ValueError("depth must be an int >= 0.")
    if isinstance(uids, str):
        uid_list = [uids]
    elif isinstance(uids, list):
        uid_list = list(uids)
    else:
        raise TypeError("uids must be a str or list[str].")
    if not uid_list:
        raise ValueError("uids is required and must contain at least one uid.")
    model_name = _registry_class_name(registry)
    queryset = getattr(DB(source), model_name)
    saved: list[SQLRecord] = []
    for uid in uid_list:
        record = queryset.get(uid)
        kwargs: dict[str, Any] = {"depth": depth}
        if transfer is not None:
            kwargs["transfer"] = transfer
        saved.append(record.save(**kwargs))
    return saved
