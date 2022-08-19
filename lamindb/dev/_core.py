from lndb_setup import settings
from lnschema_core import dobject


def storage_key_from_dobject(dobj: dobject):
    return f"{dobj.id}-{dobj.v}{dobj.file_suffix}"


def storage_key_from_triple(dobj_id: str, dobj_v: str, dobj_suffix: str):
    return f"{dobj_id}-{dobj_v}{dobj_suffix}"


def filepath_from_dobject(dobj: dobject):
    storage_key = storage_key_from_dobject(dobj)
    filepath = settings.instance.storage.key_to_filepath(storage_key)
    return filepath
