import lnschema_core as core
from lndb_setup import settings

from lamindb.dev.db._add import add
from lamindb.dev.db._select import select


def add_dobject_from_dtransform(dobject: core.dobject, dtransform_id: str):
    storage = select(core.storage, root=str(settings.instance.storage_root)).one()

    dobject_id = add(
        core.dobject(
            id=dobject.id,
            name=dobject.name,
            dtransform_id=dtransform_id,
            suffix=dobject.suffix,
            storage_id=storage.id,
            size=dobject.size,
            checksum=dobject.checksum,
        )
    )

    return dobject_id
