import sqlmodel as sqm
from lndb_setup import settings
from lnschema_core import dobject, type

from lamindb import schema


def storage_key_from_dobject(dobj: dobject):
    return f"{dobj.id}-{dobj.v}{dobj.file_suffix}"


def storage_key_from_triple(dobj_id: str, dobj_v: str, dobj_suffix: str):
    return f"{dobj_id}-{dobj_v}{dobj_suffix}"


def filepath_from_dobject(dobj: dobject):
    storage_key = storage_key_from_dobject(dobj)
    filepath = settings.instance.storage.key_to_filepath(storage_key)
    return filepath


def track_usage(dobject_id, dobject_v, usage_type: type.usage):
    with sqm.Session(settings.instance.db_engine()) as session:
        usage = schema.core.usage(
            type=usage_type,
            user_id=settings.user.id,
            dobject_id=dobject_id,
            dobject_v=dobject_v,
        )
        session.add(usage)
        session.commit()
        session.refresh(usage)

    settings.instance._update_cloud_sqlite_file()

    return usage.id
